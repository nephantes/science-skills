import type { ValidatorParams, ValidatorResult } from '../types';

/**
 * Output Validator — now a thin wrapper.
 *
 * All substantive checks (empty DataFrames, all-NaN, missing files based on
 * step title keywords) have been consolidated into `validateOutputContract()`
 * in `services/output-contract.ts` to eliminate duplicate [DIAG] parsing.
 *
 * This validator is kept as a no-op pass-through so existing skill registrations
 * don't break. The pipeline calls `validateOutputContract()` directly.
 */
export function validate(params: ValidatorParams): ValidatorResult {
  // All checks now handled by validateOutputContract() in the pipeline.
  // Returning pass to avoid duplicate issue reporting.
  void params;
  return { passed: true, issues: [] };
}
