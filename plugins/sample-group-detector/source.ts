import type { PluginParams, PluginResult } from '../types';
import * as fs from 'fs';
import * as path from 'path';

// ─── Controlled Vocabulary + Abbreviation Expansion ────────────────────────

const BIO_MAP: Record<string, Array<[string, string]>> = {
  genotype: [
    ['wt', 'Wild Type'], ['ko', 'Knockout'], ['het', 'Heterozygous'],
    ['wildtype', 'Wild Type'], ['knockout', 'Knockout'], ['heterozygous', 'Heterozygous'],
    ['mutant', 'Mutant'], ['mut', 'Mutant'], ['transgenic', 'Transgenic'],
    ['tg', 'Transgenic'], ['flox', 'Floxed'], ['cre', 'Cre Recombinase'],
    ['kd', 'Knockdown'], ['oe', 'Overexpression'], ['null', 'Null Allele'],
    ['ki', 'Knock-In'], ['wt1', 'Wild Type Rep1'], ['ko1', 'Knockout Rep1'],
    ['wt2', 'Wild Type Rep2'], ['ko2', 'Knockout Rep2'],
  ],
  diet: [
    ['chow', 'Chow Diet'], ['hfd', 'High Fat Diet'], ['hf', 'High Fat'],
    ['lfd', 'Low Fat Diet'], ['lfat', 'Low Fat'], ['hfat', 'High Fat'],
    ['normal', 'Normal Diet'], ['obesogenic', 'Obesogenic Diet'],
    ['western', 'Western Diet'], ['ketogenic', 'Ketogenic Diet'], ['keto', 'Ketogenic'],
    ['fasting', 'Fasting'], ['fed', 'Fed'], ['refed', 'Refed'],
    ['hfhs', 'High Fat High Sucrose'], ['nd', 'Normal Diet'],
    ['cd', 'Control Diet'], ['wd', 'Western Diet'],
  ],
  perturbation: [
    ['lps', 'Lipopolysaccharide'], ['ifn', 'Interferon'], ['tnf', 'TNF-alpha'],
    ['tgfb', 'TGF-beta'], ['il6', 'Interleukin-6'], ['il1b', 'Interleukin-1beta'],
    ['dex', 'Dexamethasone'], ['rapa', 'Rapamycin'], ['torin', 'Torin'],
    ['crispr', 'CRISPR'], ['sirna', 'siRNA'], ['shrna', 'shRNA'],
    ['radiation', 'Radiation'], ['irradiation', 'Irradiation'],
    ['hypoxia', 'Hypoxia'], ['normoxia', 'Normoxia'],
    ['stim', 'Stimulated'], ['unstim', 'Unstimulated'],
    ['infection', 'Infection'], ['infected', 'Infected'],
  ],
  condition: [
    ['ctrl', 'Control'], ['control', 'Control'], ['treat', 'Treated'],
    ['treated', 'Treated'], ['treatment', 'Treatment'], ['vehicle', 'Vehicle'],
    ['veh', 'Vehicle'], ['untreated', 'Untreated'], ['mock', 'Mock'],
    ['sham', 'Sham'], ['drug', 'Drug'], ['dmso', 'DMSO'],
    ['stimulated', 'Stimulated'], ['unstimulated', 'Unstimulated'],
    ['exposed', 'Exposed'], ['unexposed', 'Unexposed'],
    ['placebo', 'Placebo'], ['baseline', 'Baseline'],
  ],
  sex: [
    ['m', 'Male'], ['f', 'Female'], ['male', 'Male'], ['female', 'Female'],
  ],
  tissue: [
    ['liver', 'Liver'], ['kidney', 'Kidney'], ['brain', 'Brain'],
    ['heart', 'Heart'], ['lung', 'Lung'], ['muscle', 'Muscle'],
    ['spleen', 'Spleen'], ['adipose', 'Adipose Tissue'],
    ['cortex', 'Cortex'], ['hippocampus', 'Hippocampus'],
    ['hypothalamus', 'Hypothalamus'], ['colon', 'Colon'], ['skin', 'Skin'],
    ['blood', 'Blood'], ['serum', 'Serum'], ['plasma', 'Plasma'],
    ['pbmc', 'PBMCs'], ['bone', 'Bone'], ['marrow', 'Bone Marrow'],
    ['thymus', 'Thymus'], ['pancreas', 'Pancreas'], ['intestine', 'Intestine'],
    ['retina', 'Retina'], ['testis', 'Testis'], ['ovary', 'Ovary'],
  ],
};

/** Flat lookup: token -> [category, expansion] for fast matching. */
const TOKEN_LOOKUP = new Map<string, [string, string]>();
for (const [category, pairs] of Object.entries(BIO_MAP)) {
  for (const [token, expansion] of pairs) {
    TOKEN_LOOKUP.set(token, [category, expansion]);
  }
}

// ─── Noise Filtering ──────────────────────────────────────────────────────

const NOISE_PATTERNS = [
  /^(lane|l)[0-9]+$/i,
  /^(batch|b)[0-9]+$/i,
  /^(plate|plt?)[0-9]+$/i,
  /^[ACGT]{6,}$/,
  /^s[0-9]{3,}$/i,
  /^(lib|library)[0-9]+$/i,
  /^(run|flowcell|fc)[0-9]+$/i,
  /^[A-H](0[1-9]|1[0-2])$/,
];

function isNoiseToken(token: string): boolean {
  return NOISE_PATTERNS.some(p => p.test(token));
}

// ─── Numeric Pattern Differentiation ──────────────────────────────────────

const TIMEPOINT_PATTERNS = [
  /^(\d+)(h|hr|hrs|hour|hours)$/i,
  /^(\d+)(d|day|days)$/i,
  /^(\d+)(w|wk|wks|week|weeks)$/i,
  /^(\d+)(m|mo|mon|month|months)$/i,
  /^(\d+)(min|mins|minutes?)$/i,
  /^[tT](\d+)$/,
  /^(day|d)(\d+)$/i,
  /^(week|wk|w)(\d+)$/i,
  /^(hour|hr|h)(\d+)$/i,
];

const DOSAGE_PATTERNS = [
  /^(\d+\.?\d*)(mg|ug|ng|pg|ml|ul|mm|um|nm|pm)$/i,
  /^(\d+\.?\d*)(mgkg|mg\/kg|mpk)$/i,
  /^(low|mid|high)dose$/i,
  /^dose(\d+)$/i,
];

function classifyNumericToken(token: string): 'timepoint' | 'dosage' | 'replicate' | null {
  if (TIMEPOINT_PATTERNS.some(p => p.test(token))) return 'timepoint';
  if (DOSAGE_PATTERNS.some(p => p.test(token))) return 'dosage';
  if (isReplicateToken(token)) return 'replicate';
  return null;
}

// ─── Gene symbol patterns ─────────────────────────────────────────────────

function parseGeneModification(token: string): [string, string] | null {
  const deltaMatch = token.match(/^[\u0394\u03B4](.+)$/);
  if (deltaMatch) return [deltaMatch[1], 'deletion'];

  const suffixMatch = token.match(/^([A-Za-z0-9]+?)(KO|KD|Mut|OE|KI|Null|Het|Cre|Flox)$/i);
  if (suffixMatch && suffixMatch[1].length >= 2) {
    return [suffixMatch[1], suffixMatch[2].toUpperCase()];
  }

  return null;
}

// ─── Shared helpers ────────────────────────────────────────────────────────

const METADATA_NAMES = new Set([
  'metadata.csv', 'sample_info.csv', 'samples.csv', 'coldata.csv',
  'metadata.tsv', 'sample_info.tsv', 'samples.tsv', 'phenodata.csv',
  'sample_metadata.csv', 'sampleinfo.csv', 'clinical.csv', 'clinical.tsv',
  // .txt variants
  'metadata.txt', 'sample_info.txt', 'samples.txt', 'phenodata.txt',
  'coldata.txt', 'samplesheet.txt', 'sample_metadata.txt',
]);

function isReplicateToken(token: string): boolean {
  return /^(rep|r|Rep|R)?[0-9]+$/.test(token);
}

// ─── Detection result types ──────────────────────────────────────────────

interface DetectedGroup {
  index?: number;
  label: string;
  category: string;
  values: string[];
  confidence: 'high' | 'medium' | 'low';
  expansions?: string[];
  method: string;
}

interface FileResult {
  separator?: string;
  groups: DetectedGroup[];
  sampleCount: number;
  sampleExamples: string[];
  noiseTokens?: string[];
  factorialStructure?: string;
}

// ─── Main execute function ──────────────────────────────────────────────

// ─── Metadata file inspector ──────────────────────────────────────────────

function inspectMetadataFile(filePath: string, fileName: string, language: 'python' | 'r'): string {
  try {
    const fd = fs.openSync(filePath, 'r');
    const buf = Buffer.alloc(8192);
    const bytesRead = fs.readSync(fd, buf, 0, 8192, 0);
    fs.closeSync(fd);
    const content = buf.toString('utf-8', 0, bytesRead);
    const rawLines = content.split('\n').filter(l => l.trim().length > 0);
    if (rawLines.length < 2) {
      return `\u2713 Metadata file "${fileName}" detected \u2014 use it directly for sample grouping.`;
    }

    const delim = (rawLines[0].match(/\t/g) ?? []).length > (rawLines[0].match(/,/g) ?? []).length ? '\t' : ',';
    const sep = delim === '\t' ? '\\t' : ',';
    const headers = rawLines[0].split(delim).map(h => h.trim().replace(/^["']|["']$/g, ''));
    const dataRows = rawLines.slice(1, Math.min(rawLines.length, 35))
      .map(l => l.split(delim).map(v => v.trim().replace(/^["']|["']$/g, '')));
    const nSamples = dataRows.length;

    // Identify sample ID column (first column, or column named sample/samples/id)
    const sampleIdColNames = new Set(['sample', 'samples', 'id', 'sampleid', 'sample_id', 'samplename', 'sample_name']);
    const sampleColIdx = headers.findIndex(h => sampleIdColNames.has(h.toLowerCase()));
    const sampleColName = sampleColIdx >= 0 ? headers[sampleColIdx] : headers[0];

    interface ColInfo { name: string; uniqueNonNA: number; naCount: number; values: string[] }
    const colInfos: ColInfo[] = headers.map((h, i) => {
      const vals = dataRows.map(r => r[i] ?? '').filter(v => v !== '');
      const nonNA = vals.filter(v => v.toLowerCase() !== 'na' && v !== 'null' && v !== '');
      return {
        name: h,
        uniqueNonNA: new Set(nonNA).size,
        naCount: vals.length - nonNA.length,
        values: [...new Set(nonNA)].slice(0, 10),
      };
    });

    // Group columns: 2–15 unique non-NA values, <60% NA, not the sample ID column
    const groupCols = colInfos
      .filter(c =>
        c.name !== sampleColName &&
        c.uniqueNonNA >= 2 &&
        c.uniqueNonNA <= 15 &&
        c.naCount / nSamples <= 0.6,
      )
      // Sort: fewer NAs first, then more unique values (more granular groups preferred for DE)
      .sort((a, b) => {
        const naDiff = a.naCount - b.naCount;
        if (naDiff !== 0) return naDiff;
        return b.uniqueNonNA - a.uniqueNonNA;
      });

    if (groupCols.length === 0) {
      return `\u2713 Metadata file "${fileName}" detected (${nSamples} samples) \u2014 use it directly for sample grouping. Sample column: "${sampleColName}".`;
    }

    const lines: string[] = [
      `\u2713 Metadata file "${fileName}" detected (${nSamples} samples).`,
      `  Read with: sep="${sep}". Sample column: "${sampleColName}".`,
      `  Use this file for sample grouping \u2014 do NOT parse expression column names.`,
      `  Group columns found:`,
    ];

    for (const c of groupCols) {
      const naNote = c.naCount > 0 ? `, ${c.naCount} NA` : '';
      const vals = c.values.slice(0, 8).join(', ');
      lines.push(`    "${c.name}": ${c.uniqueNonNA} groups \u2014 ${vals}${c.values.length > 8 ? ', \u2026' : ''}${naNote}`);
    }

    if (groupCols.length >= 2) {
      // Multiple group columns — advise merging into one combined group
      const colNames = groupCols.map(c => `"${c.name}"`).join(', ');
      lines.push(`  GROUP MERGING RULE: Multiple group columns detected (${colNames}). ALWAYS merge them into a SINGLE combined "group" column for ALL analyses (DE, PCA, clustering, heatmaps, enrichment, etc.).`);
      if (language === 'r') {
        const pasteArgs = groupCols.map(c => `meta$${c.name}`).join(', ');
        lines.push(`  → meta$group <- paste(${pasteArgs}, sep="_")`);
      } else {
        const concatArgs = groupCols.map(c => `meta['${c.name}']`).join(` + '_' + `);
        lines.push(`  → meta['group'] = ${concatArgs}`);
      }
      lines.push(`  Use this merged "group" for PCA coloring, DE design/contrasts, heatmap annotations, clustering, and any analysis needing sample groups. Only keep groups separate if the user explicitly asks.`);
    } else {
      const best = groupCols[0];
      lines.push(`  \u2192 Recommended primary group column: "${best.name}" (${best.uniqueNonNA} groups: ${best.values.join(', ')})`);
    }
    lines.push(`  \u2192 In Step 1: read "${fileName}", rename "${sampleColName}" \u2192 "sample", save sample_metadata.csv.`);

    if (language === 'r') {
      lines.push(`  R read: meta <- read.csv('${fileName}', sep='${sep === '\\t' ? '\\t' : ','}', stringsAsFactors=FALSE)`);
    } else {
      lines.push(`  Python read: meta = pd.read_csv('${fileName}', sep='${sep === '\\t' ? '\\t' : ','}')`);
    }

    return lines.join('\n');
  } catch {
    return `\u2713 Metadata file "${fileName}" detected \u2014 use it directly for sample grouping.`;
  }
}

export async function execute(params: PluginParams): Promise<PluginResult> {
  // Check if a metadata file is already attached — inspect it to identify group columns
  const metaFileName = params.fileNames.find(f => METADATA_NAMES.has(f.toLowerCase()));
  if (metaFileName) {
    const metaPath = path.join(params.workDir, 'uploads', metaFileName);
    const advice = inspectMetadataFile(metaPath, metaFileName, params.language);
    return { contextText: advice };
  }

  if (params.fileNames.length === 0) {
    return { contextText: '' };
  }

  const reports: string[] = [];
  const metadata: Record<string, unknown> = {};

  for (const fileName of params.fileNames) {
    const filePath = path.join(params.workDir, 'uploads', fileName);
    if (!fs.existsSync(filePath)) continue;

    const ext = fileName.split('.').pop()?.toLowerCase() ?? '';
    if (!['csv', 'tsv', 'txt'].includes(ext)) continue;

    // Read header + a few data rows
    const fd = fs.openSync(filePath, 'r');
    const buf = Buffer.alloc(8192);
    const bytesRead = fs.readSync(fd, buf, 0, 8192, 0);
    fs.closeSync(fd);
    const head = buf.toString('utf-8', 0, bytesRead);
    const lines = head.split('\n').filter(l => l.trim().length > 0);
    if (lines.length < 2) continue;

    const delimiter = (lines[0].match(/\t/g) ?? []).length > (lines[0].match(/,/g) ?? []).length ? '\t' : ',';
    const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));

    // Skip first column (usually gene IDs/names) — rest are sample names
    const sampleNames = headers.slice(1);
    if (sampleNames.length < 2) continue;

    const result = analyzeSampleNames(sampleNames, params.language);
    if (result && result.groups.length > 0) {
      reports.push(formatReport(fileName, result, params.language));
      metadata[fileName] = result;
    }
  }

  if (reports.length === 0) {
    return { contextText: '', metadata };
  }

  // Check if multiple non-replicate groups were detected across files
  const allGroups = Object.values(metadata as Record<string, FileResult>)
    .flatMap(r => r.groups)
    .filter(g => g.category !== 'unknown' || g.values.length >= 2);
  const namedGroupCount = allGroups.filter(g =>
    g.category !== 'unknown' && g.category !== 'timepoint' && g.category !== 'dosage',
  ).length;

  const mergeAdvice = namedGroupCount >= 2
    ? [
        '',
        `GROUP MERGING RULE: Multiple group variables were detected. ALWAYS merge them into a SINGLE combined "group" column using paste() (R) or string concatenation (Python).`,
        `For example, if genotype and diet are detected, create: group = paste(genotype, diet, sep="_") → "wt_chow", "ko_hfd", etc.`,
        `Use this single merged "group" column for ALL downstream analyses: PCA coloring, DE design/contrasts, heatmap annotations, clustering labels, enrichment, and any plot needing sample groups.`,
        `Only keep groups separate if the user explicitly asks to test individual factors or interactions.`,
      ]
    : [];

  const contextText = [
    `\uD83D\uDD0D SAMPLE GROUP DETECTION (no metadata file found):`,
    ...reports,
    '',
    `IMPORTANT: In Step 1, you MUST parse these groups from the column names and save them to sample_metadata.csv.`,
    `Do NOT hardcode group assignments \u2014 derive them programmatically from the sample names using the separator and positions shown above.`,
    `Include ALL detected variables as columns in sample_metadata.csv.`,
    ...mergeAdvice,
  ].join('\n');

  return { contextText, metadata };
}

// ─── Core analysis engine ──────────────────────────────────────────────

function analyzeSampleNames(sampleNames: string[], _language: 'python' | 'r'): FileResult | null {
  const separators = ['.', '_', '-'];
  let bestSep = '';
  let bestScore = 0;
  for (const sep of separators) {
    const score = sampleNames.filter(s => s.includes(sep)).length / sampleNames.length;
    if (score > bestScore) {
      bestScore = score;
      bestSep = sep;
    }
  }

  if (bestScore < 0.5 || !bestSep) {
    return analyzeWholeNames(sampleNames);
  }

  const parts = sampleNames.map(s => s.split(bestSep));
  const maxParts = Math.max(...parts.map(p => p.length));
  const partMatrix = parts.map(p => {
    const padded = [...p];
    while (padded.length < maxParts) padded.push('');
    return padded;
  });

  const detectedGroups: DetectedGroup[] = [];
  const usedIndices = new Set<number>();
  const noiseTokens: string[] = [];

  // Noise Filtering
  for (let col = 0; col < maxParts; col++) {
    const colValues = partMatrix.map(row => row[col]).filter(v => v !== '');
    if (colValues.length > 0 && colValues.every(isNoiseToken)) {
      noiseTokens.push(...[...new Set(colValues)]);
      usedIndices.add(col);
    }
  }

  // Pass 1: Controlled Vocabulary Matching
  // Two thresholds: ≥50% match → high/medium confidence (strict),
  // ≥1 match with remaining tokens as bio identifiers → medium/low (lenient).
  // Lenient match catches e.g. genotype positions with gene names like [dbl, J1c, J2c, wt]
  // where only "wt" is in the vocabulary but the position clearly represents genotype.
  for (let col = 0; col < maxParts; col++) {
    if (usedIndices.has(col)) continue;
    const colValues = partMatrix.map(row => row[col].toLowerCase());
    const uniqueValues = [...new Set(colValues)].filter(v => v !== '');
    if (uniqueValues.length === 0) continue;

    for (const [category, pairs] of Object.entries(BIO_MAP)) {
      const knownTokens = pairs.map(p => p[0]);
      const matchCount = uniqueValues.filter(v => knownTokens.includes(v)).length;
      if (matchCount > 0 && matchCount >= uniqueValues.length * 0.5) {
        // Strict match: ≥50% tokens recognized
        const expansions = uniqueValues.map(v => {
          const entry = TOKEN_LOOKUP.get(v);
          return entry ? entry[1] : v;
        });
        detectedGroups.push({
          index: col, label: category, category, values: uniqueValues,
          confidence: matchCount === uniqueValues.length ? 'high' : 'medium',
          expansions, method: 'vocabulary',
        });
        usedIndices.add(col);
        break;
      } else if (
        matchCount > 0 &&
        uniqueValues.length >= 2 &&
        uniqueValues.length <= 10 &&
        // Remaining unmatched tokens must look like biological identifiers
        // (alphanumeric, 2+ chars, not pure replicates or noise)
        uniqueValues
          .filter(v => !knownTokens.includes(v))
          .every(v => /^[a-zA-Z][a-zA-Z0-9]{1,}$/.test(v) && !isReplicateToken(v) && !isNoiseToken(v))
      ) {
        // Lenient match: at least one known token + remaining look like bio IDs (e.g. gene names)
        const expansions = uniqueValues.map(v => {
          const entry = TOKEN_LOOKUP.get(v);
          return entry ? entry[1] : v;
        });
        detectedGroups.push({
          index: col, label: category, category, values: uniqueValues,
          confidence: 'low',
          expansions, method: 'vocabulary-partial',
        });
        usedIndices.add(col);
        break;
      }
    }
  }

  // Pass 2: Timepoint + Numeric Pattern Recognition
  for (let col = 0; col < maxParts; col++) {
    if (usedIndices.has(col)) continue;
    const colValues = partMatrix.map(row => row[col]).filter(v => v !== '');
    const uniqueValues = [...new Set(colValues)];
    if (uniqueValues.length === 0) continue;

    const classifications = uniqueValues.map(classifyNumericToken);
    const timepointCount = classifications.filter(c => c === 'timepoint').length;
    const dosageCount = classifications.filter(c => c === 'dosage').length;
    const replicateCount = classifications.filter(c => c === 'replicate').length;

    if (timepointCount > 0 && timepointCount >= uniqueValues.length * 0.5) {
      detectedGroups.push({
        index: col, label: 'timepoint', category: 'timepoint', values: uniqueValues,
        confidence: timepointCount === uniqueValues.length ? 'high' : 'medium',
        method: 'pattern',
      });
      usedIndices.add(col);
    } else if (dosageCount > 0 && dosageCount >= uniqueValues.length * 0.5) {
      detectedGroups.push({
        index: col, label: 'dosage', category: 'dosage', values: uniqueValues,
        confidence: dosageCount === uniqueValues.length ? 'high' : 'medium',
        method: 'pattern',
      });
      usedIndices.add(col);
    } else if (replicateCount === uniqueValues.length) {
      usedIndices.add(col);
    }
  }

  // Pass 3: Gene Symbol + Modification Patterns
  for (let col = 0; col < maxParts; col++) {
    if (usedIndices.has(col)) continue;
    const colValues = partMatrix.map(row => row[col]).filter(v => v !== '');
    const uniqueValues = [...new Set(colValues)];
    if (uniqueValues.length === 0) continue;

    const geneModifications = uniqueValues.map(parseGeneModification).filter(Boolean);
    if (geneModifications.length > 0 && geneModifications.length >= uniqueValues.length * 0.5) {
      const expansions = uniqueValues.map(v => {
        const parsed = parseGeneModification(v);
        return parsed ? `${parsed[0]}-${parsed[1]}` : v;
      });
      detectedGroups.push({
        index: col, label: 'genotype', category: 'genotype', values: uniqueValues,
        confidence: 'medium', expansions, method: 'gene-symbol',
      });
      usedIndices.add(col);
    }
  }

  // Pass 4: Differential Token Analysis + Frequency Clustering
  for (let col = 0; col < maxParts; col++) {
    if (usedIndices.has(col)) continue;
    const colValues = partMatrix.map(row => row[col]).filter(v => v !== '');
    const uniqueValues = [...new Set(colValues)];

    if (uniqueValues.length <= 1) continue;
    if (uniqueValues.length >= sampleNames.length * 0.9) continue;
    if (uniqueValues.every(isReplicateToken)) continue;
    if (uniqueValues.every(isNoiseToken)) continue;

    const freq = new Map<string, number>();
    for (const v of colValues) {
      freq.set(v, (freq.get(v) ?? 0) + 1);
    }
    const freqValues = [...freq.values()];
    const isBalanced = freqValues.length >= 2 &&
      Math.max(...freqValues) / Math.min(...freqValues) <= 3;

    let confidence: 'high' | 'medium' | 'low' = 'low';
    if (isBalanced && uniqueValues.length >= 2 && uniqueValues.length <= 6) {
      confidence = 'medium';
    }

    const looksLikeCondition = uniqueValues.some(v =>
      /^[a-zA-Z]{2,}[0-9]*$/.test(v) && !isReplicateToken(v),
    );
    if (looksLikeCondition) {
      confidence = confidence === 'low' ? 'medium' : confidence;
    }

    detectedGroups.push({
      index: col, label: `group_${col + 1}`, category: 'unknown',
      values: uniqueValues, confidence, method: 'frequency',
    });
  }

  if (detectedGroups.length === 0) return null;

  const factorialStructure = validateFactorialStructure(detectedGroups, partMatrix);

  return {
    separator: bestSep,
    groups: detectedGroups,
    sampleCount: sampleNames.length,
    sampleExamples: sampleNames.slice(0, 4),
    noiseTokens: noiseTokens.length > 0 ? noiseTokens : undefined,
    factorialStructure,
  };
}

// ─── Whole-name analysis (no separator) ──────────────────────────────────

function analyzeWholeNames(sampleNames: string[]): FileResult | null {
  const lower = sampleNames.map(s => s.toLowerCase());
  const groups: DetectedGroup[] = [];

  for (const [category, pairs] of Object.entries(BIO_MAP)) {
    const knownTokens = pairs.map(p => p[0]);
    const matchingValues = [...new Set(lower)].filter(v =>
      knownTokens.some(t => v.includes(t)),
    );
    if (matchingValues.length >= 2) {
      const expansions = matchingValues.map(v => {
        for (const [token, expansion] of pairs) {
          if (v.includes(token)) return expansion;
        }
        return v;
      });
      groups.push({
        label: category, category, values: matchingValues,
        confidence: 'medium', expansions, method: 'vocabulary-substring',
      });
    }
  }

  const timepointMatches = [...new Set(lower)].filter(v =>
    TIMEPOINT_PATTERNS.some(p => p.test(v)),
  );
  if (timepointMatches.length >= 2) {
    groups.push({
      label: 'timepoint', category: 'timepoint', values: timepointMatches,
      confidence: 'high', method: 'pattern',
    });
  }

  const geneMatches = sampleNames.filter(s => parseGeneModification(s) !== null);
  if (geneMatches.length >= 2) {
    const uniqueGene = [...new Set(geneMatches)];
    groups.push({
      label: 'genotype', category: 'genotype', values: uniqueGene,
      confidence: 'medium',
      expansions: uniqueGene.map(v => {
        const parsed = parseGeneModification(v);
        return parsed ? `${parsed[0]}-${parsed[1]}` : v;
      }),
      method: 'gene-symbol',
    });
  }

  if (groups.length === 0) return null;
  return {
    groups,
    sampleCount: sampleNames.length,
    sampleExamples: sampleNames.slice(0, 4),
  };
}

// ─── Factorial Structure Validation ─────────────────────────────────────

function validateFactorialStructure(
  groups: DetectedGroup[],
  partMatrix: string[][],
): string | undefined {
  const indexedGroups = groups.filter(g => g.index !== undefined);
  if (indexedGroups.length < 2) return undefined;

  const combinations = new Set<string>();
  for (const row of partMatrix) {
    const combo = indexedGroups.map(g => row[g.index!]).join(' \u00D7 ');
    combinations.add(combo);
  }

  let expectedCount = 1;
  for (const g of indexedGroups) {
    expectedCount *= g.values.length;
  }

  const observedCount = combinations.size;
  const factorLabels = indexedGroups.map(g => `${g.label}(${g.values.length})`).join(' \u00D7 ');

  if (observedCount === expectedCount) {
    return `Full factorial design: ${factorLabels} = ${expectedCount} combinations (all present)`;
  } else if (observedCount >= expectedCount * 0.7) {
    return `Near-complete factorial: ${factorLabels} \u2014 ${observedCount}/${expectedCount} combinations observed`;
  } else {
    return `Partial factorial: ${factorLabels} \u2014 ${observedCount}/${expectedCount} combinations (some missing)`;
  }
}

// ─── Report formatting ──────────────────────────────────────────────────

function formatReport(
  fileName: string,
  result: FileResult,
  language: 'python' | 'r',
): string {
  const lines: string[] = [];
  lines.push(`\nFile "${fileName}" \u2014 ${result.sampleCount} samples (e.g. ${result.sampleExamples.join(', ')})`);

  if (result.separator) {
    lines.push(`  Separator: "${result.separator}"`);
  }

  if (result.noiseTokens && result.noiseTokens.length > 0) {
    lines.push(`  Noise tokens (ignored): ${result.noiseTokens.join(', ')}`);
  }

  for (const g of result.groups) {
    const posInfo = g.index !== undefined ? ` [position ${g.index + 1}]` : '';
    const confIcon = g.confidence === 'high' ? '\u2713' : g.confidence === 'medium' ? '~' : '?';
    const methodInfo = g.method !== 'vocabulary' ? ` (via ${g.method})` : '';

    if (g.expansions && g.expansions.length > 0) {
      const expanded = g.values.map((v, i) => `${v} (${g.expansions![i]})`).join(', ');
      lines.push(`  ${confIcon} ${g.label}: ${expanded}${posInfo}${methodInfo}`);
    } else {
      lines.push(`  ${confIcon} ${g.label}: ${g.values.join(', ')}${posInfo}${methodInfo}`);
    }
  }

  if (result.factorialStructure) {
    lines.push(`  Design: ${result.factorialStructure}`);
  }

  if (result.separator) {
    const indexedGroups = result.groups.filter(g => g.index !== undefined);
    const groupPositions = indexedGroups.map(g => g.index! + 1); // 1-based

    const slotMap = indexedGroups
      .map(g => `position ${g.index! + 1} \u2192 ${g.label}`)
      .join(', ');
    if (slotMap) {
      lines.push(`  Slot mapping: ${slotMap} (replicate suffix excluded from group)`);
    }

    // Generate concrete extraction code — explicit enough that Claude cannot misread the separator
    if (groupPositions.length > 0) {
      const exampleGroup = result.sampleExamples[0]
        .split(result.separator)
        .filter((_, i) => groupPositions.includes(i + 1))
        .join(result.separator);

      if (language === 'r') {
        const posStr = groupPositions.length === 1
          ? String(groupPositions[0])
          : `c(${groupPositions.join(', ')})`;
        lines.push(`  Group extraction (R) — USE THIS EXACT PATTERN:`);
        lines.push(`    parts <- strsplit(sample_names, "${result.separator}", fixed=TRUE)`);
        lines.push(`    group <- sapply(parts, function(x) paste(x[${posStr}], collapse="${result.separator}"))`);
        lines.push(`    # "${result.sampleExamples[0]}" → "${exampleGroup}"`);
      } else {
        const sep = result.separator;
        const pyIndices = groupPositions.map(p => p - 1); // 0-based
        const posStr = pyIndices.length === 1 ? String(pyIndices[0]) : `[${pyIndices.join(', ')}]`;
        const sliceExpr = pyIndices.length === 1
          ? `parts[${posStr}]`
          : `[parts[i] for i in ${posStr}]`;
        lines.push(`  Group extraction (Python) \u2014 USE THIS EXACT PATTERN:`);
        lines.push(`    parts_list = [c.split('${sep}') for c in df.columns]`);
        lines.push(`    group = ['${sep}'.join(${sliceExpr.replace('parts', 'p')}) for p in parts_list]`);
        lines.push(`    # "${result.sampleExamples[0]}" \u2192 "${exampleGroup}"`);
      }
    }
  }

  return lines.join('\n');
}
