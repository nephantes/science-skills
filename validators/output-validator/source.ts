import type { ValidatorParams, ValidatorResult } from '../types';

/** Step title keywords that imply file output is expected. */
const FILE_PRODUCING_KEYWORDS = [
  'plot', 'save', 'figure', 'export', 'visualization',
  'heatmap', 'volcano', 'chart', 'graph', 'write', 'generate',
];

/**
 * Parse [DIAG] lines from stdout to extract DataFrame inventory.
 * Format: [DIAG] df 'name': ROWSxCOLS, NaN=PCT%, cols=[...]
 */
function parseDiagDataFrames(stdout: string): Array<{
  name: string;
  rows: number;
  cols: number;
  nanPct: number;
}> {
  const dfs: Array<{ name: string; rows: number; cols: number; nanPct: number }> = [];
  const regex = /\[DIAG\] df '(\w+)': (\d+)x(\d+), NaN=([\d.]+)%/g;
  let match;
  while ((match = regex.exec(stdout)) !== null) {
    dfs.push({
      name: match[1],
      rows: parseInt(match[2], 10),
      cols: parseInt(match[3], 10),
      nanPct: parseFloat(match[4]),
    });
  }
  return dfs;
}

export function validate(params: ValidatorParams): ValidatorResult {
  const issues: string[] = [];

  // Only validate successful steps — hard failures go straight to the fix loop
  if (params.exitCode !== 0) {
    return { passed: true, issues: [] };
  }

  // Check 1: Empty DataFrames (0 rows)
  const dfs = parseDiagDataFrames(params.stdout);
  for (const df of dfs) {
    if (df.rows === 0) {
      issues.push(
        `DataFrame '${df.name}' has 0 rows — filtering or merge likely removed all data. Check filter conditions and merge keys.`,
      );
    }
  }

  // Check 2: All-NaN DataFrames
  for (const df of dfs) {
    if (df.rows > 0 && df.nanPct >= 100) {
      issues.push(
        `DataFrame '${df.name}' is 100% NaN (${df.rows}x${df.cols}). Data type coercion likely failed — check if columns contain non-numeric strings.`,
      );
    }
  }

  // Check 3: No output files when title implies file production
  const titleLower = params.stepTitle.toLowerCase();
  const matchedKeywords = FILE_PRODUCING_KEYWORDS.filter(kw => titleLower.includes(kw));
  if (matchedKeywords.length > 0 && params.outputFileNames.length === 0) {
    issues.push(
      `Step "${params.stepTitle}" produced no output files but title suggests it should (matched: ${matchedKeywords.join(', ')}). Ensure plots are saved with savefig()/ggsave() and data with to_csv()/write.csv().`,
    );
  }

  return { passed: issues.length === 0, issues };
}
