import type { PluginParams, PluginResult } from '../types';

const PYTHON_RULES = `
CRITICAL CODE QUALITY RULES — SILENT ERROR PREVENTION:
1. NEVER use bare try/except or try/except Exception that silently passes:
   BAD:  try: ... except: pass
   BAD:  try: ... except Exception: continue
   GOOD: try: ... except KeyError as e: print(f"Warning: {e}", flush=True); raise
2. NEVER use warnings.filterwarnings('ignore') globally — only suppress specific known warnings.
3. ALWAYS add shape/type assertions after data transformations:
   - After merge: assert len(merged) > 0, "Merge produced empty result"
   - After filter: assert len(filtered) > 0, f"Filter removed all rows"
   - After read_csv: print(f"Loaded {df.shape[0]} rows x {df.shape[1]} cols", flush=True)
4. ALWAYS flush print output: print(..., flush=True)
5. For resource safety:
   - Limit DataFrame operations to avoid memory blowup: sample large datasets before correlation/distance matrices
   - When computing pairwise distances on >10000 features, subset first
6. ALWAYS print intermediate results (shapes, column names, group counts) so errors are diagnosable.
7. NEVER catch and silently convert exceptions to default values without printing a warning.`.trim();

const R_RULES = `
CRITICAL CODE QUALITY RULES — SILENT ERROR PREVENTION:
1. NEVER wrap entire script in tryCatch — only wrap specific risky operations.
2. NEVER use suppressWarnings() or suppressMessages() globally.
   GOOD: suppressWarnings(single_expression)  # for a known benign warning
   BAD:  suppressWarnings({ ... entire block ... })
3. ALWAYS add validation after data transformations:
   - After merge: stopifnot(nrow(merged) > 0)
   - After filter: if(nrow(filtered) == 0) stop("Filter removed all rows")
   - After read.csv: cat("Loaded:", nrow(df), "rows x", ncol(df), "cols\\n")
4. For resource safety:
   - Limit matrix operations: if ncol(mat) > 10000, sample columns first
   - Check memory before large operations: cat("Matrix size:", object.size(mat), "bytes\\n")
5. ALWAYS print intermediate results (dimensions, column names, factor levels).
6. NEVER use silent=TRUE in download.file or similar functions — always check return values.`.trim();

export async function execute(params: PluginParams): Promise<PluginResult> {
  const contextText = params.language === 'python' ? PYTHON_RULES : R_RULES;
  return { contextText };
}
