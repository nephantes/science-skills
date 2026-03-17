// Re-export skill types so plugin source.ts files can import from '../types'
// This barrel resolves the relative import path from mnt/skills-repo/plugins/*/source.ts
export type { PluginParams, PluginResult } from '../../../src/skills/types';
