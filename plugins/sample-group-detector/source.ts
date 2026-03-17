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
  timepoint: [
    // Timepoint tokens are detected via regex, not vocabulary lookup.
    // This entry exists so labeling can assign "timepoint" as a category.
  ],
  dosage: [
    // Dosage tokens are detected via regex, not vocabulary lookup.
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

function isTimepointToken(token: string): boolean {
  return TIMEPOINT_PATTERNS.some(p => p.test(token));
}

function isDosageToken(token: string): boolean {
  return DOSAGE_PATTERNS.some(p => p.test(token));
}

// ─── Replicate detection ─────────────────────────────────────────────────

/**
 * Check if a single token matches a replicate pattern.
 * Matches: rep1, r1, R1, Rep1, 1, 01, 001, etc.
 */
function isReplicateToken(token: string): boolean {
  return /^(rep|r|Rep|R)\d+$/.test(token) || /^\d+$/.test(token);
}

/**
 * Check if ALL values in a column match replicate patterns.
 * This is the strict dataset-level check.
 */
function isReplicateColumn(values: string[]): boolean {
  if (values.length === 0) return false;
  return values.every(v => isReplicateToken(v));
}

/**
 * Check if values look sequential/numeric (characteristic of replicates).
 * Returns true if values are sequential integers or rep-prefixed sequential integers.
 */
function areValuesSequentialOrNumeric(values: string[]): boolean {
  // Extract numeric parts
  const nums = values.map(v => {
    const m = v.match(/^(?:rep|r|Rep|R)?(\d+)$/);
    return m ? parseInt(m[1], 10) : NaN;
  });
  if (nums.some(isNaN)) return false;
  // Check if they form a reasonable range (not all the same, max - min < count * 3)
  const unique = [...new Set(nums)];
  if (unique.length <= 1) return false;
  const min = Math.min(...unique);
  const max = Math.max(...unique);
  return max - min < unique.length * 3;
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

// ─── Detection result types ──────────────────────────────────────────────

interface ColumnLabel {
  position: number;
  category: string;
  label: string;
  values: string[];
  expansion?: string[];
}

interface GroupingCandidate {
  /** Column indices that form the group (0-based) */
  groupColumns: number[];
  /** Column indices identified as replicates (0-based) */
  replicateColumns: number[];
  /** Column indices identified as noise (0-based) */
  noiseColumns: number[];
  /** Inferred group labels per sample */
  groups: string[];
  /** Unique group labels */
  uniqueGroups: string[];
  /** Replicate values per sample */
  replicates: string[];
  /** Confidence score */
  confidence: 'high' | 'medium' | 'low';
  /** Scoring breakdown */
  score: number;
}

interface FileResult {
  separator: string;
  grouping: GroupingCandidate;
  columnLabels: ColumnLabel[];
  sampleCount: number;
  sampleExamples: string[];
  noiseTokens?: string[];
  factorialStructure?: string;
  timepointDetected?: { position: number; values: string[] };
}

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

// ─── Main execute function ──────────────────────────────────────────────

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
    if (result && result.grouping.uniqueGroups.length > 0) {
      reports.push(formatReport(fileName, result, params.language));
      metadata[fileName] = result;
    }
  }

  if (reports.length === 0) {
    return { contextText: '', metadata };
  }

  // Check if multiple non-replicate group variables were detected
  const allResults = Object.values(metadata) as FileResult[];
  const allLabels = allResults.flatMap(r => r.columnLabels);
  const namedGroupCount = allLabels.filter(l =>
    l.category !== 'replicate' && l.category !== 'noise' && l.category !== 'unknown',
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

// ─── Core analysis engine (dataset-level inference) ────────────────────────

function analyzeSampleNames(sampleNames: string[], _language: 'python' | 'r'): FileResult | null {
  // Step 1: Detect dominant separator
  const separator = detectDominantSeparator(sampleNames);
  if (!separator) {
    return analyzeWholeNames(sampleNames);
  }

  // Step 2: Tokenize all names, build token matrix
  const tokenized = sampleNames.map(s => s.split(separator));
  const maxTokens = Math.max(...tokenized.map(t => t.length));

  // Pad shorter arrays to build a uniform matrix
  const matrix: string[][] = tokenized.map(tokens => {
    const padded = [...tokens];
    while (padded.length < maxTokens) padded.push('');
    return padded;
  });

  // Step 3: Identify replicate columns (dataset-level)
  const replicateColumns = detectReplicateColumns(matrix, sampleNames.length);

  // Step 4: Identify noise columns
  const noiseColumns = detectNoiseColumns(matrix);

  // Step 5: Remaining columns form the group
  const allColumns = Array.from({ length: maxTokens }, (_, i) => i);
  const groupColumns = allColumns.filter(
    c => !replicateColumns.includes(c) && !noiseColumns.includes(c),
  );

  // Build group labels by joining non-replicate, non-noise tokens
  const groups = matrix.map(row =>
    groupColumns.map(c => row[c]).filter(v => v !== '').join(separator),
  );
  const uniqueGroups = [...new Set(groups)];
  const replicates = matrix.map(row =>
    replicateColumns.map(c => row[c]).filter(v => v !== '').join(separator),
  );

  // Step 6: Candidate scoring
  const candidate = scoreCandidate(
    groupColumns, replicateColumns, noiseColumns,
    groups, uniqueGroups, replicates, sampleNames.length,
  );

  // If ambiguous, try alternative candidates
  const bestCandidate = resolveAmbiguity(candidate, matrix, separator, sampleNames.length);

  // Step 7: Vocabulary labeling (cosmetic enrichment)
  const columnLabels = labelColumns(matrix, separator);

  // Detect timepoint columns
  let timepointDetected: { position: number; values: string[] } | undefined;
  for (const label of columnLabels) {
    if (label.category === 'timepoint') {
      timepointDetected = { position: label.position, values: label.values };
      break;
    }
  }

  // Collect noise tokens
  const noiseTokens: string[] = [];
  for (const col of bestCandidate.noiseColumns) {
    const colVals = matrix.map(row => row[col]).filter(v => v !== '');
    noiseTokens.push(...[...new Set(colVals)]);
  }

  // Factorial structure
  const factorialStructure = computeFactorialStructure(
    columnLabels.filter(l =>
      bestCandidate.groupColumns.includes(l.position) &&
      l.category !== 'unknown' && l.category !== 'noise',
    ),
    matrix,
  );

  return {
    separator,
    grouping: bestCandidate,
    columnLabels,
    sampleCount: sampleNames.length,
    sampleExamples: sampleNames.slice(0, 4),
    noiseTokens: noiseTokens.length > 0 ? noiseTokens : undefined,
    factorialStructure,
    timepointDetected,
  };
}

// ─── Step 1: Separator detection ──────────────────────────────────────────

function detectDominantSeparator(sampleNames: string[]): string | null {
  const separators = ['.', '_', '-'];
  const counts: Record<string, number> = {};

  for (const sep of separators) {
    counts[sep] = 0;
    for (const name of sampleNames) {
      // Count how many names contain this separator
      if (name.includes(sep)) counts[sep]++;
    }
  }

  // Find the separator used in the most names
  let bestSep = '';
  let bestCount = 0;
  for (const sep of separators) {
    if (counts[sep] > bestCount) {
      bestCount = counts[sep];
      bestSep = sep;
    }
  }

  // Require at least 50% of samples to have the separator
  const ratio = bestCount / sampleNames.length;
  if (ratio < 0.5 || !bestSep) return null;

  // Handle mixed separators: if a secondary separator is used in >20% of names,
  // normalize those names by replacing the secondary with the dominant separator
  // (This is handled during tokenization, not here — we just return the dominant one.)

  return bestSep;
}

// ─── Step 3: Replicate column detection (dataset-level) ───────────────────

function detectReplicateColumns(matrix: string[][], sampleCount: number): number[] {
  const numCols = matrix[0].length;
  const candidates: Array<{ col: number; score: number }> = [];

  for (let col = 0; col < numCols; col++) {
    const colValues = matrix.map(row => row[col]).filter(v => v !== '');
    if (colValues.length === 0) continue;

    // Check if ALL values match replicate patterns
    if (!isReplicateColumn(colValues)) continue;

    const unique = [...new Set(colValues)];
    const cardinality = unique.length;

    // Cardinality-to-sample-count ratio: replicates have high cardinality relative to
    // the number of times each value repeats (each value appears few times = many groups use it)
    // BUT also replicates typically have LOW cardinality relative to total samples
    // (e.g., 3 replicates across 6 samples = cardinality 3, each appears 2 times)

    // Key heuristic: replicate column has values that are sequential/numeric
    const sequential = areValuesSequentialOrNumeric(unique);

    // Score: prefer columns that look like replicates
    let score = 0;
    if (sequential) score += 10;
    // Replicates typically have low-to-medium cardinality (2-10)
    if (cardinality >= 2 && cardinality <= 10) score += 5;
    // Each replicate value should appear roughly the same number of times
    const freq = new Map<string, number>();
    for (const v of colValues) freq.set(v, (freq.get(v) ?? 0) + 1);
    const freqs = [...freq.values()];
    const maxFreq = Math.max(...freqs);
    const minFreq = Math.min(...freqs);
    if (maxFreq > 0 && minFreq > 0 && maxFreq / minFreq <= 2) score += 3;
    // Prefer last positions (replicates are usually at the end)
    score += col * 0.5;

    if (score > 0) {
      candidates.push({ col, score });
    }
  }

  if (candidates.length === 0) return [];

  // If only one candidate, use it
  if (candidates.length === 1) return [candidates[0].col];

  // If multiple candidates, pick the one with the highest score
  // Typically there should be only one replicate column
  candidates.sort((a, b) => b.score - a.score);
  return [candidates[0].col];
}

// ─── Step 4: Noise column detection ───────────────────────────────────────

function detectNoiseColumns(matrix: string[][]): number[] {
  const numCols = matrix[0].length;
  const noiseCols: number[] = [];

  for (let col = 0; col < numCols; col++) {
    const colValues = matrix.map(row => row[col]).filter(v => v !== '');
    if (colValues.length > 0 && colValues.every(isNoiseToken)) {
      noiseCols.push(col);
    }
  }

  return noiseCols;
}

// ─── Step 6: Candidate scoring ────────────────────────────────────────────

function scoreCandidate(
  groupColumns: number[],
  replicateColumns: number[],
  noiseColumns: number[],
  groups: string[],
  uniqueGroups: string[],
  replicates: string[],
  sampleCount: number,
): GroupingCandidate {
  let score = 0;
  let confidence: 'high' | 'medium' | 'low' = 'low';

  // 1. No singletons (every group should appear more than once)
  const groupFreq = new Map<string, number>();
  for (const g of groups) groupFreq.set(g, (groupFreq.get(g) ?? 0) + 1);
  const freqs = [...groupFreq.values()];
  const singletonCount = freqs.filter(f => f === 1).length;
  const hasSingletons = singletonCount > 0;

  if (!hasSingletons) {
    score += 20;
  } else if (singletonCount <= uniqueGroups.length * 0.3) {
    score += 10; // A few singletons are okay (incomplete design)
  }

  // 2. Balanced replicate counts
  if (freqs.length > 0) {
    const maxFreq = Math.max(...freqs);
    const minFreq = Math.min(...freqs);
    if (maxFreq > 0 && minFreq > 0) {
      const ratio = maxFreq / minFreq;
      if (ratio <= 1.5) score += 15;
      else if (ratio <= 2.0) score += 10;
      else if (ratio <= 3.0) score += 5;
    }
  }

  // 3. Repeated groups (groups should appear multiple times if replicates exist)
  if (replicateColumns.length > 0 && uniqueGroups.length < sampleCount) {
    score += 10;
  }

  // 4. Reasonable number of groups
  if (uniqueGroups.length >= 2 && uniqueGroups.length <= Math.ceil(sampleCount / 2)) {
    score += 10;
  }

  // 5. Group columns exist
  if (groupColumns.length > 0) {
    score += 5;
  }

  // 6. Replicate columns detected
  if (replicateColumns.length > 0) {
    score += 5;
  }

  // Determine confidence
  if (score >= 45 && !hasSingletons && replicateColumns.length > 0) {
    confidence = 'high';
  } else if (score >= 25) {
    confidence = 'medium';
  } else {
    confidence = 'low';
  }

  return {
    groupColumns,
    replicateColumns,
    noiseColumns,
    groups,
    uniqueGroups,
    replicates,
    confidence,
    score,
  };
}

// ─── Ambiguity resolution ─────────────────────────────────────────────────

/**
 * If the primary candidate is ambiguous, try alternative interpretations
 * and return the best one.
 */
function resolveAmbiguity(
  primary: GroupingCandidate,
  matrix: string[][],
  separator: string,
  sampleCount: number,
): GroupingCandidate {
  // If primary is high confidence, no need to explore alternatives
  if (primary.confidence === 'high') return primary;

  const numCols = matrix[0].length;

  // If no replicates were detected, try treating each column as a potential replicate
  if (primary.replicateColumns.length === 0) {
    let bestAlt: GroupingCandidate | null = null;

    for (let col = 0; col < numCols; col++) {
      if (primary.noiseColumns.includes(col)) continue;
      const colValues = matrix.map(row => row[col]).filter(v => v !== '');
      if (colValues.length === 0) continue;

      // Try this column as replicate even if not all values match strict patterns
      // (e.g., some datasets use arbitrary sample IDs)
      const unique = [...new Set(colValues)];
      if (unique.length < 2) continue;
      if (unique.length >= sampleCount * 0.8) continue; // Too many unique = likely replicate

      // Build alternative grouping
      const altGroupCols = Array.from({ length: numCols }, (_, i) => i).filter(
        c => c !== col && !primary.noiseColumns.includes(c),
      );
      const altGroups = matrix.map(row =>
        altGroupCols.map(c => row[c]).filter(v => v !== '').join(separator),
      );
      const altUniqueGroups = [...new Set(altGroups)];
      const altReplicates = matrix.map(row => row[col]);

      const altCandidate = scoreCandidate(
        altGroupCols, [col], primary.noiseColumns,
        altGroups, altUniqueGroups, altReplicates, sampleCount,
      );

      if (!bestAlt || altCandidate.score > bestAlt.score) {
        bestAlt = altCandidate;
      }
    }

    if (bestAlt && bestAlt.score > primary.score) {
      return bestAlt;
    }
  }

  return primary;
}

// ─── Step 7: Vocabulary labeling (cosmetic enrichment) ────────────────────

function labelColumns(matrix: string[][], _separator: string): ColumnLabel[] {
  const numCols = matrix[0].length;
  const labels: ColumnLabel[] = [];

  for (let col = 0; col < numCols; col++) {
    const colValues = matrix.map(row => row[col]).filter(v => v !== '');
    const unique = [...new Set(colValues)];
    if (unique.length === 0) continue;

    // Check noise first
    if (unique.every(isNoiseToken)) {
      labels.push({ position: col, category: 'noise', label: 'noise', values: unique });
      continue;
    }

    // Check if replicate column
    if (isReplicateColumn(colValues) && areValuesSequentialOrNumeric(unique)) {
      labels.push({ position: col, category: 'replicate', label: 'replicate', values: unique });
      continue;
    }

    // Check timepoint patterns
    if (unique.some(isTimepointToken) && unique.filter(isTimepointToken).length >= unique.length * 0.5) {
      labels.push({ position: col, category: 'timepoint', label: 'timepoint', values: unique });
      continue;
    }

    // Check dosage patterns
    if (unique.some(isDosageToken) && unique.filter(isDosageToken).length >= unique.length * 0.5) {
      labels.push({ position: col, category: 'dosage', label: 'dosage', values: unique });
      continue;
    }

    // Check gene modification patterns
    const geneMatches = unique.filter(v => parseGeneModification(v) !== null);
    if (geneMatches.length > 0 && geneMatches.length >= unique.length * 0.5) {
      const expansion = unique.map(v => {
        const parsed = parseGeneModification(v);
        return parsed ? `${parsed[0]}-${parsed[1]}` : v;
      });
      labels.push({ position: col, category: 'genotype', label: 'genotype', values: unique, expansion });
      continue;
    }

    // Vocabulary matching (BIO_MAP)
    let matched = false;
    const lowerUnique = unique.map(v => v.toLowerCase());
    for (const [category, pairs] of Object.entries(BIO_MAP)) {
      if (pairs.length === 0) continue; // skip empty entries (timepoint, dosage)
      const knownTokens = pairs.map(p => p[0]);
      const matchCount = lowerUnique.filter(v => knownTokens.includes(v)).length;

      if (matchCount > 0 && matchCount >= lowerUnique.length * 0.5) {
        // Strict match
        const expansion = lowerUnique.map(v => {
          const entry = TOKEN_LOOKUP.get(v);
          return entry ? entry[1] : v;
        });
        labels.push({ position: col, category, label: category, values: unique, expansion });
        matched = true;
        break;
      } else if (
        matchCount > 0 &&
        lowerUnique.length >= 2 &&
        lowerUnique.length <= 10 &&
        lowerUnique
          .filter(v => !knownTokens.includes(v))
          .every(v => /^[a-zA-Z][a-zA-Z0-9]{1,}$/.test(v) && !isReplicateToken(v) && !isNoiseToken(v))
      ) {
        // Lenient match: at least one known token + remaining look biological
        const expansion = lowerUnique.map(v => {
          const entry = TOKEN_LOOKUP.get(v);
          return entry ? entry[1] : v;
        });
        labels.push({ position: col, category, label: category, values: unique, expansion });
        matched = true;
        break;
      }
    }

    if (!matched) {
      // Unknown column — label as generic group
      labels.push({ position: col, category: 'unknown', label: `position_${col + 1}`, values: unique });
    }
  }

  return labels;
}

// ─── Factorial structure ──────────────────────────────────────────────────

function computeFactorialStructure(
  groupLabels: ColumnLabel[],
  matrix: string[][],
): string | undefined {
  if (groupLabels.length < 2) return undefined;

  const combinations = new Set<string>();
  for (const row of matrix) {
    const combo = groupLabels.map(g => row[g.position]).join(' \u00D7 ');
    combinations.add(combo);
  }

  let expectedCount = 1;
  for (const g of groupLabels) {
    expectedCount *= g.values.length;
  }

  const observedCount = combinations.size;
  const factorLabels = groupLabels.map(g => `${g.label}(${g.values.length})`).join(' \u00D7 ');

  if (observedCount === expectedCount) {
    return `Full factorial design: ${factorLabels} = ${expectedCount} combinations (all present)`;
  } else if (observedCount >= expectedCount * 0.7) {
    return `Near-complete factorial: ${factorLabels} \u2014 ${observedCount}/${expectedCount} combinations observed`;
  } else {
    return `Partial factorial: ${factorLabels} \u2014 ${observedCount}/${expectedCount} combinations (some missing)`;
  }
}

// ─── Whole-name analysis (no separator) ──────────────────────────────────

function analyzeWholeNames(sampleNames: string[]): FileResult | null {
  // Try to find a common prefix/suffix pattern
  // e.g., control1, control2, treated1, treated2
  const patterns = detectPrefixSuffixGroups(sampleNames);
  if (patterns) return patterns;

  // Fall back to vocabulary substring matching
  const lower = sampleNames.map(s => s.toLowerCase());
  const columnLabels: ColumnLabel[] = [];

  for (const [category, pairs] of Object.entries(BIO_MAP)) {
    if (pairs.length === 0) continue;
    const knownTokens = pairs.map(p => p[0]);
    const matchingValues = [...new Set(lower)].filter(v =>
      knownTokens.some(t => v.includes(t)),
    );
    if (matchingValues.length >= 2) {
      const expansion = matchingValues.map(v => {
        for (const [token, exp] of pairs) {
          if (v.includes(token)) return exp;
        }
        return v;
      });
      columnLabels.push({
        position: 0, label: category, category, values: matchingValues, expansion,
      });
    }
  }

  // Gene modification patterns
  const geneMatches = sampleNames.filter(s => parseGeneModification(s) !== null);
  if (geneMatches.length >= 2) {
    const uniqueGene = [...new Set(geneMatches)];
    columnLabels.push({
      position: 0, label: 'genotype', category: 'genotype', values: uniqueGene,
      expansion: uniqueGene.map(v => {
        const parsed = parseGeneModification(v);
        return parsed ? `${parsed[0]}-${parsed[1]}` : v;
      }),
    });
  }

  if (columnLabels.length === 0) return null;

  // Build a basic grouping from vocabulary matches
  const groups = sampleNames.map(s => {
    // Strip trailing digits to get group
    const m = s.match(/^(.*?)(\d+)$/);
    return m ? m[1] : s;
  });
  const uniqueGroups = [...new Set(groups)];

  return {
    separator: '',
    grouping: {
      groupColumns: [0],
      replicateColumns: [],
      noiseColumns: [],
      groups,
      uniqueGroups,
      replicates: sampleNames.map(s => {
        const m = s.match(/(\d+)$/);
        return m ? m[1] : '';
      }),
      confidence: 'low',
      score: 10,
    },
    columnLabels,
    sampleCount: sampleNames.length,
    sampleExamples: sampleNames.slice(0, 4),
  };
}

/**
 * Detect groups from names without separators by finding a common structure
 * where a prefix varies and a numeric suffix acts as replicate.
 * e.g., control1, control2, treated1, treated2
 */
function detectPrefixSuffixGroups(sampleNames: string[]): FileResult | null {
  // Try splitting each name into alpha prefix + numeric suffix
  const parsed = sampleNames.map(s => {
    const m = s.match(/^(.*?)(\d+)$/);
    if (!m) return null;
    return { prefix: m[1], num: m[2], original: s };
  });

  // If not all names have this pattern, bail
  if (parsed.some(p => p === null)) return null;
  const valid = parsed as Array<{ prefix: string; num: string; original: string }>;

  const prefixes = [...new Set(valid.map(p => p.prefix))];
  const nums = [...new Set(valid.map(p => p.num))];

  // Need at least 2 samples and some structure
  if (prefixes.length < 1 || nums.length < 1) return null;
  // If every name has a unique prefix, no grouping is possible
  if (prefixes.length === valid.length) return null;

  const groups = valid.map(p => p.prefix);
  const uniqueGroups = [...new Set(groups)];
  const replicates = valid.map(p => p.num);

  const grouping = scoreCandidate(
    [0], [1], [],
    groups, uniqueGroups, replicates, sampleNames.length,
  );

  // Label the prefix using vocabulary
  const columnLabels: ColumnLabel[] = [];
  const lowerPrefixes = prefixes.map(p => p.toLowerCase());
  let labeled = false;
  for (const [category, pairs] of Object.entries(BIO_MAP)) {
    if (pairs.length === 0) continue;
    const knownTokens = pairs.map(p => p[0]);
    const matchCount = lowerPrefixes.filter(v => knownTokens.includes(v)).length;
    if (matchCount > 0) {
      const expansion = lowerPrefixes.map(v => {
        const entry = TOKEN_LOOKUP.get(v);
        return entry ? entry[1] : v;
      });
      columnLabels.push({ position: 0, category, label: category, values: prefixes, expansion });
      labeled = true;
      break;
    }
  }
  if (!labeled) {
    columnLabels.push({ position: 0, category: 'unknown', label: 'group', values: prefixes });
  }
  columnLabels.push({ position: 1, category: 'replicate', label: 'replicate', values: nums });

  return {
    separator: '',
    grouping,
    columnLabels,
    sampleCount: sampleNames.length,
    sampleExamples: sampleNames.slice(0, 4),
  };
}

// ─── Report formatting ──────────────────────────────────────────────────

function formatReport(
  fileName: string,
  result: FileResult,
  language: 'python' | 'r',
): string {
  const lines: string[] = [];
  const { grouping, columnLabels, separator } = result;

  lines.push(`\nFile "${fileName}" \u2014 ${result.sampleCount} samples (e.g. ${result.sampleExamples.join(', ')})`);

  if (separator) {
    lines.push(`  Separator: "${separator}"`);
  }

  if (result.noiseTokens && result.noiseTokens.length > 0) {
    lines.push(`  Noise tokens (ignored): ${result.noiseTokens.join(', ')}`);
  }

  // Show detected group structure
  lines.push(`  Groups detected: ${grouping.uniqueGroups.join(', ')} (${grouping.uniqueGroups.length} groups)`);
  lines.push(`  Confidence: ${grouping.confidence}`);

  // Show replicate info
  if (grouping.replicateColumns.length > 0) {
    const repValues = [...new Set(grouping.replicates)].filter(v => v !== '');
    lines.push(`  Replicate pattern: ${repValues.join(', ')} [position ${grouping.replicateColumns.map(c => c + 1).join(', ')}]`);
  }

  // Show column labels (vocabulary enrichment)
  const meaningfulLabels = columnLabels.filter(l =>
    l.category !== 'replicate' && l.category !== 'noise' && l.category !== 'unknown',
  );
  for (const label of meaningfulLabels) {
    const posInfo = ` [position ${label.position + 1}]`;
    if (label.expansion && label.expansion.length > 0) {
      const expanded = label.values.map((v, i) => `${v} (${label.expansion![i]})`).join(', ');
      lines.push(`  \u2713 ${label.label}: ${expanded}${posInfo}`);
    } else {
      lines.push(`  \u2713 ${label.label}: ${label.values.join(', ')}${posInfo}`);
    }
  }

  // Show example mapping
  if (result.sampleExamples.length > 0) {
    const example = result.sampleExamples[0];
    const exampleIdx = 0;
    const exampleGroup = grouping.groups[exampleIdx];
    const exampleRep = grouping.replicates[exampleIdx];
    lines.push(`  Example: "${example}" \u2192 group="${exampleGroup}", replicate="${exampleRep}"`);
  }

  if (result.timepointDetected) {
    lines.push(`  Timepoint column detected [position ${result.timepointDetected.position + 1}]: ${result.timepointDetected.values.join(', ')}`);
  }

  if (result.factorialStructure) {
    lines.push(`  Design: ${result.factorialStructure}`);
  }

  // Generic extraction instructions (language-agnostic)
  if (separator) {
    const groupPositions = grouping.groupColumns.map(c => c + 1); // 1-based
    const repPosition = grouping.replicateColumns.length > 0 ? grouping.replicateColumns[0] + 1 : undefined;

    const slotMap = columnLabels
      .filter(l => l.category !== 'noise')
      .map(l => `position ${l.position + 1} = ${l.label}`)
      .join(', ');
    if (slotMap) {
      lines.push(`  Parsing rule: split each sample name on "${separator}" (literal character, not regex) → ${slotMap}`);
    }

    lines.push(`  Group = join tokens at position${groupPositions.length > 1 ? 's' : ''} ${groupPositions.join(', ')} with "${separator}"`);
    if (repPosition) {
      lines.push(`  Replicate = token at position ${repPosition}`);
    }

    // Show individual variable columns for multi-factor designs
    const namedLabels = columnLabels.filter(l =>
      l.category !== 'replicate' && l.category !== 'noise' && l.category !== 'unknown',
    );
    if (namedLabels.length >= 2) {
      lines.push(`  Save each factor as a separate column in sample_metadata.csv:`);
      for (const label of namedLabels) {
        lines.push(`    ${label.label} = token at position ${label.position + 1}`);
      }
    }
  } else {
    lines.push(`  Parsing rule: strip trailing digits to get group, extract trailing digits as replicate`);
  }

  return lines.join('\n');
}
