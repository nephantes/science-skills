import type { PluginParams, PluginResult } from '../types';
import * as fs from 'fs';
import * as path from 'path';

// ─── Column name patterns (lowercased for matching) ───────────────────────────

/** Log fold change column names */
const LFC_PATTERNS = new Set([
  'log2foldchange', 'log2fc', 'logfc', 'lfc', 'log_fc',
  'log2_fold_change', 'log2fold', 'log2_fc', 'foldchange', 'fold_change',
  'logfoldchange', 'log.fold.change', 'log2.fold.change',
]);

/** Adjusted p-value column names */
const PADJ_PATTERNS = new Set([
  'padj', 'adj.p.val', 'adj_p_val', 'adj.pval', 'adj_pval',
  'fdr', 'q.value', 'qvalue', 'q_value', 'p.adjust', 'p_adjust',
  'adjusted_pvalue', 'adjusted.pvalue', 'adjustedpvalue',
  'bh', 'fdr_bh', 'p.adj', 'p_adj',
]);

/** Comparison / contrast column names */
const COMPARISON_PATTERNS = new Set([
  'comparison', 'contrast', 'group_comparison', 'group_contrast',
  'comparison_name', 'contrast_name', 'de_comparison', 'test',
]);

// ─── Helpers ──────────────────────────────────────────────────────────────────

interface DeFileInfo {
  fileName: string;
  lfcCol: string;
  padjCol: string;
  /** Name of a column containing comparison labels, if present */
  comparisonCol: string | null;
  /** Unique values in the comparison column (up to 20) */
  comparisons: string[];
  /** Prefixed comparison blocks found in wide format (e.g. "KO_vs_WT") */
  prefixedComparisons: string[];
  nRows: number;
}

/** Detect prefixed multi-comparison columns (wide format).
 *  e.g. log2FoldChange_KO_vs_WT and padj_KO_vs_WT → comparison "KO_vs_WT"
 */
function detectPrefixedComparisons(headers: string[], origHeaders: string[]): string[] {
  const comparisons = new Set<string>();

  for (const h of headers) {
    for (const lfcPat of LFC_PATTERNS) {
      if (h.startsWith(lfcPat + '_') || h.startsWith(lfcPat + '.')) {
        const suffix = h.slice(lfcPat.length + 1);
        if (suffix.length > 0) {
          // Verify there's a matching padj column with the same suffix
          const hasPadj = headers.some(oh =>
            PADJ_PATTERNS.has(oh.replace(new RegExp(`[_.]${suffix.replace(/[-/]/g, '[-/]')}$`), ''))
          );
          if (hasPadj || headers.some(oh => oh.includes('padj') && oh.includes(suffix))) {
            // Recover original casing from origHeaders
            const origIdx = headers.indexOf(h);
            const orig = origHeaders[origIdx] ?? h;
            const origSuffix = orig.slice(orig.indexOf('_') + 1);
            comparisons.add(origSuffix || suffix);
          }
        }
      }
    }
  }

  return [...comparisons];
}

/** Read the first N unique values of a column from parsed rows */
function readColumnValues(
  rows: string[][],
  colIdx: number,
  maxValues: number,
): string[] {
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

    // Read up to 32 KB to parse header + sample rows
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

    // Find LFC and padj columns (exact match first)
    const lfcIdx = headers.findIndex(h => LFC_PATTERNS.has(h));
    const padjIdx = headers.findIndex(h => PADJ_PATTERNS.has(h));

    if (lfcIdx === -1 || padjIdx === -1) continue;

    const lfcCol = origHeaders[lfcIdx];
    const padjCol = origHeaders[padjIdx];
    const nRows = rawLines.length - 1;

    // Check for a comparison/contrast column
    const compIdx = headers.findIndex(h => COMPARISON_PATTERNS.has(h));
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
      // Long format with a comparison column
      lines.push(`  Comparison column  : "${info.comparisonCol}" — ${info.comparisons.length} comparison(s):`);
      for (const c of info.comparisons) {
        lines.push(`    • ${c}`);
      }
      lines.push(`  Format: LONG (filter by "${info.comparisonCol}" to get per-comparison results)`);
    } else if (info.prefixedComparisons.length > 1) {
      // Wide format with prefixed columns
      lines.push(`  Wide-format comparisons detected (${info.prefixedComparisons.length}):`);
      for (const c of info.prefixedComparisons) {
        lines.push(`    • ${c}`);
      }
      lines.push(`  Format: WIDE (each comparison has its own log2FC/padj column pair)`);
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
      lines.push(`   Loop over each column pair (log2FC_<comparison>, padj_<comparison>) to generate per-comparison plots.`);
    }
    lines.push(`   Also produce a summary plot (DEG counts barplot or heatmap) across all comparisons.`);
  } else {
    lines.push(`3. Produce a volcano plot (x = ${detected[0].lfcCol}, y = -log10(${detected[0].padjCol})) and any other requested downstream analyses.`);
  }

  return { contextText: lines.join('\n') };
}
