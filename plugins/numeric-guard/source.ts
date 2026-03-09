import type { PluginParams, PluginResult } from '../types';
import * as fs from 'fs';
import * as path from 'path';

/** Keywords in user queries that indicate numeric operations are needed. */
const NUMERIC_KEYWORDS = [
  'pca', 'clustering', 'heatmap', 'correlation', 'normalization', 'normalize',
  'differential expression', 'deseq2', 'statistics', 'regression', 'anova',
  'ttest', 't-test', 'fold change', 'volcano', 'scatter', 'boxplot',
  'variance', 'standard deviation', 'mean', 'median', 'distribution',
  'dimensionality reduction', 'umap', 'tsne', 'enrichment', 'pathway',
  'gene expression', 'rna-seq', 'rnaseq', 'counts', 'fpkm', 'tpm',
  'survival', 'kaplan', 'cox', 'distance', 'similarity',
];

/** Detect whether a value is likely numeric (int or float). */
function isNumericValue(v: string): boolean {
  const trimmed = v.trim();
  if (trimmed === '' || trimmed === 'NA' || trimmed === 'NaN' || trimmed === 'null') return false;
  return !isNaN(Number(trimmed));
}

export async function execute(params: PluginParams): Promise<PluginResult> {
  const queryLower = params.userMessage.toLowerCase();
  const needsNumeric = NUMERIC_KEYWORDS.some(kw => queryLower.includes(kw));

  const fileReports: string[] = [];
  const metadata: Record<string, unknown> = {};

  for (const fileName of params.fileNames) {
    const filePath = path.join(params.workDir, 'uploads', fileName);
    if (!fs.existsSync(filePath)) continue;

    const ext = fileName.split('.').pop()?.toLowerCase() ?? '';
    if (!['csv', 'tsv', 'txt'].includes(ext)) continue;

    // Read header + first data rows
    const fd = fs.openSync(filePath, 'r');
    const buf = Buffer.alloc(8192);
    const bytesRead = fs.readSync(fd, buf, 0, 8192, 0);
    fs.closeSync(fd);
    const head = buf.toString('utf-8', 0, bytesRead);
    const lines = head.split('\n').filter(l => l.trim().length > 0);
    if (lines.length < 2) continue;

    const delimiter = (lines[0].match(/\t/g) ?? []).length > (lines[0].match(/,/g) ?? []).length ? '\t' : ',';
    const headers = lines[0].split(delimiter).map(h => h.trim().replace(/^["']|["']$/g, ''));

    // Sample up to 5 data rows to classify each column
    const dataRows = lines.slice(1, 6);
    const numericCols: string[] = [];
    const nonNumericCols: string[] = [];

    for (let col = 0; col < headers.length; col++) {
      const values = dataRows.map(row => {
        const cells = row.split(delimiter);
        return (cells[col] ?? '').trim().replace(/^["']|["']$/g, '');
      }).filter(v => v !== '' && v !== 'NA' && v !== 'NaN');

      if (values.length === 0) {
        numericCols.push(headers[col]);
      } else if (values.every(isNumericValue)) {
        numericCols.push(headers[col]);
      } else {
        nonNumericCols.push(headers[col]);
      }
    }

    metadata[fileName] = {
      delimiter: delimiter === '\t' ? 'tab' : 'comma',
      totalColumns: headers.length,
      numericCols,
      nonNumericCols,
    };

    if (nonNumericCols.length > 0) {
      fileReports.push(
        `File "${fileName}" has ${nonNumericCols.length} non-numeric column(s): ${nonNumericCols.join(', ')}\n` +
        `  Numeric columns (${numericCols.length}): ${numericCols.slice(0, 10).join(', ')}${numericCols.length > 10 ? ` ... (+${numericCols.length - 10} more)` : ''}`,
      );
    }
  }

  if (fileReports.length === 0) {
    return { contextText: '', metadata };
  }

  const dropInstructions = params.language === 'r'
    ? `In R: after reading the data, select only numeric columns:\n` +
      `  df <- df[, sapply(df, is.numeric)]\n` +
      `For DESeq2/edgeR: ensure the count matrix contains ONLY integer count columns — drop any text/annotation columns first.`
    : `In Python: after reading the data, select only numeric columns:\n` +
      `  df = df.select_dtypes(include='number')\n` +
      `For statistical analysis: ensure the matrix contains ONLY numeric columns — drop any text/annotation columns first.`;

  const contextText = [
    `\u26a0\ufe0f NON-NUMERIC COLUMN GUARD:`,
    ...fileReports,
    '',
    needsNumeric
      ? `This analysis requires numeric operations. You MUST drop non-numeric columns before any computation (PCA, DESeq2, clustering, correlation, etc.).`
      : `If your analysis involves numeric operations, drop non-numeric columns before computation.`,
    dropInstructions,
  ].join('\n');

  return { contextText, metadata };
}
