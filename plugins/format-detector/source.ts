import type { PluginParams, PluginResult, ValidatorParams, ValidatorResult } from '../types';
import * as fs from 'fs';
import * as path from 'path';

export async function execute(params: PluginParams): Promise<PluginResult> {
  const detections: string[] = [];
  const metadata: Record<string, unknown> = {};

  for (const fileName of params.fileNames) {
    const filePath = path.join(params.workDir, 'uploads', fileName);
    if (!fs.existsSync(filePath)) continue;

    const ext = fileName.split('.').pop()?.toLowerCase() ?? '';
    const stat = fs.statSync(filePath);
    const sizeKB = (stat.size / 1024).toFixed(1);

    // Read first few lines for heuristics
    const fd = fs.openSync(filePath, 'r');
    const buf = Buffer.alloc(4096);
    const bytesRead = fs.readSync(fd, buf, 0, 4096, 0);
    fs.closeSync(fd);
    const head = buf.toString('utf-8', 0, bytesRead);
    const lines = head.split('\n').slice(0, 5);
    const firstLine = lines[0] ?? '';

    // Detect delimiter
    const tabCount = (firstLine.match(/\t/g) ?? []).length;
    const commaCount = (firstLine.match(/,/g) ?? []).length;
    const delimiter = tabCount > commaCount ? 'tab' : 'comma';
    const columnCount = firstLine.split(delimiter === 'tab' ? '\t' : ',').length;

    // Detect if it looks like a count matrix (integer-heavy)
    const dataLine = lines[1] ?? '';
    const values = dataLine.split(delimiter === 'tab' ? '\t' : ',').slice(1);
    const isCountMatrix = values.length > 5 &&
      values.slice(0, 10).every(v => /^\d+$/.test(v.trim()));

    detections.push(
      `- ${fileName}: ${ext.toUpperCase()}, ${sizeKB} KB, ${delimiter}-delimited, ` +
      `${columnCount} columns${isCountMatrix ? ', appears to be a count matrix' : ''}`,
    );

    metadata[fileName] = { ext, delimiter, columnCount, isCountMatrix };
  }

  return {
    contextText: detections.length > 0
      ? `File format detection results:\n${detections.join('\n')}`
      : '',
    metadata,
  };
}

// ─── Validator ──────────────────────────────────────────────────────────────

const CSV_TSV_EXTS = new Set(['csv', 'tsv', 'txt']);

export function validate(params: ValidatorParams): ValidatorResult {
  const issues: string[] = [];

  // Only validate successful steps
  if (params.exitCode !== 0) {
    return { passed: true, issues: [] };
  }

  for (const fileName of params.outputFileNames) {
    const ext = fileName.split('.').pop()?.toLowerCase() ?? '';
    if (!CSV_TSV_EXTS.has(ext)) continue;

    const filePath = path.join(params.workDir, 'outputs', fileName);
    if (!fs.existsSync(filePath)) continue;

    // Read first 8 KB for heuristics
    const fd = fs.openSync(filePath, 'r');
    const buf = Buffer.alloc(8192);
    const bytesRead = fs.readSync(fd, buf, 0, 8192, 0);
    fs.closeSync(fd);

    if (bytesRead === 0) {
      issues.push(`Output file "${fileName}" is empty (0 bytes).`);
      continue;
    }

    const head = buf.toString('utf-8', 0, bytesRead);
    const lines = head.split('\n').filter((l) => l.trim().length > 0);

    // Check 1: Must have at least a header + 1 data row
    if (lines.length < 2) {
      issues.push(
        `Output file "${fileName}" has only ${lines.length} non-empty line(s) — expected at least a header row and one data row.`,
      );
      continue;
    }

    // Detect delimiter
    const headerLine = lines[0];
    const tabCount = (headerLine.match(/\t/g) ?? []).length;
    const commaCount = (headerLine.match(/,/g) ?? []).length;
    const delimiter = tabCount > commaCount ? '\t' : ',';
    const headerCols = headerLine.split(delimiter).length;

    // Check 2: Delimiter consistency — all rows should have the same column count
    const inconsistentRows: number[] = [];
    for (let i = 1; i < Math.min(lines.length, 20); i++) {
      const cols = lines[i].split(delimiter).length;
      if (cols !== headerCols) {
        inconsistentRows.push(i + 1); // 1-indexed for humans
      }
    }
    if (inconsistentRows.length > 0) {
      issues.push(
        `Output file "${fileName}" has inconsistent column counts: header has ${headerCols} columns but row(s) ${inconsistentRows.join(', ')} differ. Check delimiter usage and quoting.`,
      );
    }

    // Check 3: Header row shouldn't be purely numeric (likely missing header)
    const headerValues = headerLine.split(delimiter);
    const allNumeric = headerValues.every((v) => /^-?\d+(\.\d+)?([eE][+-]?\d+)?$/.test(v.trim()));
    if (allNumeric && headerValues.length > 1) {
      issues.push(
        `Output file "${fileName}" appears to have no header row — all values in the first row are numeric. Ensure the code writes column names.`,
      );
    }
  }

  return { passed: issues.length === 0, issues };
}
