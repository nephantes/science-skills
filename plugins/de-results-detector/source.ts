import type { PluginParams, PluginResult } from '../types';
import * as fs from 'fs';
import * as path from 'path';

// ─── Column name matchers (all matching is case-insensitive on lowercased header) ──

/**
 * Test whether a lowercased column name is a fold-change column.
 * Covers: log2FoldChange (DESeq2), logFC/FC (edgeR/limma), lfc, fold_change,
 * foldChange, fc, L2FC, and prefixed variants like log2FC_KO_vs_WT.
 */
function isLfcColumn(h: string): boolean {
  // Exact short-form matches
  if (h === 'lfc' || h === 'l2fc') return true;
  // "fc" alone — safe because we require a padj column too
  if (h === 'fc') return true;
  // "fc" as a prefix/suffix with a separator (FC_treated, result_FC)
  if (/^fc[_.\-]/.test(h) || /[_.\-]fc$/.test(h)) return true;

  // log2 variants: log2FoldChange, log2fc, log2_fc, log2.fc, log2fold …
  if (/log2[_.\-]?fold[_.\-]?change/.test(h)) return true;
  if (/log2[_.\-]?fc/.test(h)) return true;
  if (/log2fold/.test(h)) return true;

  // logFC variants (natural log or unspecified base)
  if (/log[_.\-]?fc/.test(h)) return true;
  if (/log[_.\-]?fold[_.\-]?change/.test(h)) return true;
  if (/logfoldchange/.test(h)) return true;

  // plain fold change (linear — still a pre-computed DE result)
  if (/^fold[_.\-]?change/.test(h)) return true;
  if (/[_.\-]fold[_.\-]?change$/.test(h)) return true;
  if (/^foldchange/.test(h)) return true;

  return false;
}

/**
 * Test whether a lowercased column name is an adjusted p-value column.
 * Covers: padj (DESeq2), adj.P.Val (limma), FDR (edgeR), q.value, p.adjust, BH …
 */
function isPadjColumn(h: string): boolean {
  // Common exact matches
  if (h === 'padj' || h === 'p.adj' || h === 'p_adj') return true;
  if (h === 'fdr') return true;
  if (h === 'bh' || h === 'fdr_bh' || h === 'fdr.bh') return true;
  if (h === 'qvalue' || h === 'q.value' || h === 'q_value') return true;
  if (h === 'p.adjust' || h === 'p_adjust' || h === 'padjust') return true;

  // adj.P.Val (limma) and variants
  if (/adj[_.\-]?p[_.\-]?val/.test(h)) return true;
  // adjusted_pvalue etc.
  if (/adjusted[_.\-]?p/.test(h)) return true;
  // q.val / q_val
  if (/^q[_.\-]?val/.test(h)) return true;
  // FDR as prefix/suffix (FDR_KO_vs_WT, result_FDR)
  if (/^fdr[_.\-]/.test(h) || /[_.\-]fdr$/.test(h)) return true;
  // padj as prefix/suffix (padj_BH, result_padj)
  if (/^padj[_.\-]/.test(h) || /[_.\-]padj$/.test(h)) return true;

  return false;
}

/**
 * Returns true if the column name already implies a log-scale fold change
 * (i.e. it contains "log" or "ln" as a word component).
 */
function isLogScaleByName(h: string): boolean {
  return /log/.test(h) || /\bln\b/.test(h) || h === 'lfc' || h === 'l2fc';
}

/** Comparison / contrast column names (exact, lowercased) */
const COMPARISON_EXACT = new Set([
  'comparison', 'contrast', 'group_comparison', 'group_contrast',
  'comparison_name', 'contrast_name', 'de_comparison', 'test',
  'condition_pair', 'group_pair',
]);

// ─── Types ────────────────────────────────────────────────────────────────────

interface DeFileInfo {
  fileName: string;
  lfcCol: string;
  padjCol: string;
  comparisonCol: string | null;
  comparisons: string[];
  prefixedComparisons: string[];
  nRows: number;
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

/**
 * Detect wide-format multi-comparison columns.
 * e.g. "log2FoldChange_KO_vs_WT" + "padj_KO_vs_WT" → comparison "KO_vs_WT"
 */
function detectPrefixedComparisons(headers: string[], origHeaders: string[]): string[] {
  const comparisons = new Set<string>();

  for (let i = 0; i < headers.length; i++) {
    const h = headers[i];
    for (const sep of ['_', '.', '-']) {
      const sepIdx = h.indexOf(sep);
      if (sepIdx <= 0) continue;
      const prefix = h.slice(0, sepIdx);
      const suffix = h.slice(sepIdx + 1);
      if (!suffix) continue;

      if (isLfcColumn(prefix)) {
        // Check there's a matching padj column with the same suffix
        const hasPadjMatch = headers.some((other, j) => {
          if (j === i) return false;
          const otherSepIdx = other.indexOf(sep);
          if (otherSepIdx <= 0) return false;
          const otherPrefix = other.slice(0, otherSepIdx);
          const otherSuffix = other.slice(otherSepIdx + 1);
          return isPadjColumn(otherPrefix) && otherSuffix === suffix;
        });
        if (hasPadjMatch) {
          const orig = origHeaders[i] ?? h;
          const origSepIdx = orig.indexOf(sep);
          const origSuffix = origSepIdx >= 0 ? orig.slice(origSepIdx + 1) : suffix;
          comparisons.add(origSuffix);
        }
      }
    }
  }

  return [...comparisons];
}

function readColumnValues(rows: string[][], colIdx: number, maxValues: number): string[] {
  const seen = new Set<string>();
  for (const row of rows) {
    const val = (row[colIdx] ?? '').trim().replace(/^["']|["']$/g, '');
    if (val && val.toLowerCase() !== 'na' && val.toLowerCase() !== 'null') {
      seen.add(val);
      if (seen.size >= maxValues) break;
    }
  }
  return [...seen];
}

// ─── Main execute ──────────────────────────────────────────────────────────────

export async function execute(params: PluginParams): Promise<PluginResult> {
  if (params.fileNames.length === 0) {
    return { contextText: '' };
  }

  const detected: DeFileInfo[] = [];

  for (const fileName of params.fileNames) {
    const filePath = path.join(params.workDir, 'uploads', fileName);
    if (!fs.existsSync(filePath)) continue;

    const ext = fileName.split('.').pop()?.toLowerCase() ?? '';
    if (!['csv', 'tsv', 'txt'].includes(ext)) continue;

    // Read up to 32 KB (header + sample rows)
    const fd = fs.openSync(filePath, 'r');
    const buf = Buffer.alloc(32_768);
    const bytesRead = fs.readSync(fd, buf, 0, 32_768, 0);
    fs.closeSync(fd);

    const content = buf.toString('utf-8', 0, bytesRead);
    const rawLines = content.split('\n').filter(l => l.trim().length > 0);
    if (rawLines.length < 2) continue;

    const firstLine = rawLines[0];
    const tabCount = (firstLine.match(/\t/g) ?? []).length;
    const commaCount = (firstLine.match(/,/g) ?? []).length;
    const delimiter = tabCount > commaCount ? '\t' : ',';

    const origHeaders = firstLine.split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));
    const headers = origHeaders.map(h => h.toLowerCase());

    // Require BOTH an LFC and a padj column — either exact or prefixed
    const lfcIdx = headers.findIndex(h => isLfcColumn(h));
    const padjIdx = headers.findIndex(h => isPadjColumn(h));

    if (lfcIdx === -1 || padjIdx === -1) continue;

    const lfcCol = origHeaders[lfcIdx];
    const padjCol = origHeaders[padjIdx];
    const nRows = rawLines.length - 1;

    // Check for a comparison/contrast column (long format)
    const compIdx = headers.findIndex(h => COMPARISON_EXACT.has(h));
    let comparisonCol: string | null = null;
    let comparisons: string[] = [];

    if (compIdx !== -1) {
      comparisonCol = origHeaders[compIdx];
      const dataRows = rawLines.slice(1, Math.min(rawLines.length, 500))
        .map(l => l.split(delimiter));
      comparisons = readColumnValues(dataRows, compIdx, 20);
    }

    // Check for wide-format prefixed comparisons
    const prefixedComparisons = detectPrefixedComparisons(headers, origHeaders);

    detected.push({ fileName, lfcCol, padjCol, comparisonCol, comparisons, prefixedComparisons, nRows });
  }

  if (detected.length === 0) {
    return { contextText: '' };
  }

  // ─── Build context text ────────────────────────────────────────────────────

  const lines: string[] = [
    `⚠️  DE RESULTS DETECTOR — Pre-computed differential expression data found:`,
    ``,
  ];

  for (const info of detected) {
    lines.push(`FILE: "${info.fileName}" (${info.nRows} rows)`);
    lines.push(`  Fold-change column : "${info.lfcCol}"`);
    lines.push(`  Adjusted p-value   : "${info.padjCol}"`);

    if (info.comparisons.length > 0) {
      lines.push(`  Comparison column  : "${info.comparisonCol}" — ${info.comparisons.length} comparison(s):`);
      for (const c of info.comparisons) {
        lines.push(`    • ${c}`);
      }
      lines.push(`  Format: LONG (filter by "${info.comparisonCol}" to get per-comparison results)`);
    } else if (info.prefixedComparisons.length > 1) {
      lines.push(`  Wide-format comparisons detected (${info.prefixedComparisons.length}):`);
      for (const c of info.prefixedComparisons) {
        lines.push(`    • ${c}`);
      }
      lines.push(`  Format: WIDE (each comparison has its own fold-change/padj column pair)`);
    } else {
      lines.push(`  Format: SINGLE comparison`);
    }
    lines.push('');
  }

  lines.push(`CRITICAL PLANNING RULES:`);
  lines.push(`1. DO NOT run DESeq2, edgeR, limma, pydeseq2, or any other DE tool — results already exist.`);
  lines.push(`2. Step 1 must load the DE results file and use "${detected[0].lfcCol}" / "${detected[0].padjCol}" directly.`);

  const multiComp = detected.find(d => d.comparisons.length > 1 || d.prefixedComparisons.length > 1);
  if (multiComp) {
    const nComps = multiComp.comparisons.length > 1 ? multiComp.comparisons.length : multiComp.prefixedComparisons.length;
    lines.push(`3. MULTIPLE COMPARISONS FOUND (${nComps}): Create one volcano plot (or MA / scatter plot) per comparison.`);
    if (multiComp.comparisons.length > 1) {
      lines.push(`   Use a loop over unique values of "${multiComp.comparisonCol}" to generate per-comparison plots.`);
      lines.push(`   Save each plot with the comparison name in the filename (e.g. volcano_KO_vs_WT.png).`);
    } else {
      lines.push(`   Loop over each column pair (fold-change_<comparison>, padj_<comparison>) to generate per-comparison plots.`);
    }
    lines.push(`   Also produce a summary plot (DEG counts barplot or heatmap) across all comparisons.`);
  } else {
    lines.push(`3. Produce a volcano plot (x = ${detected[0].lfcCol}, y = -log10(${detected[0].padjCol})) and any other requested downstream analyses.`);
  }

  return { contextText: lines.join('\n') };
}
