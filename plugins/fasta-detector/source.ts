import type { PluginParams, PluginResult } from '../types';
import * as fs from 'fs';
import * as path from 'path';

const FASTA_EXTENSIONS = new Set(['fasta', 'fa', 'faa', 'fna', 'fas', 'mafft', 'aln', 'aligned']);

/** Check if a filename (possibly multi-extension like foo.faa.mafft) is FASTA-related */
function isFastaFile(fileName: string): boolean {
  const parts = fileName.split('.');
  return parts.some(p => FASTA_EXTENSIONS.has(p.toLowerCase()));
}

/** Check if content looks like FASTA (starts with ">") */
function looksLikeFasta(head: string): boolean {
  const firstNonEmpty = head.split('\n').find(l => l.trim().length > 0);
  return firstNonEmpty?.startsWith('>') ?? false;
}

export async function execute(params: PluginParams): Promise<PluginResult> {
  const fastaFiles: string[] = [];
  const detections: string[] = [];
  const metadata: Record<string, unknown> = {};

  for (const fileName of params.fileNames) {
    const filePath = path.join(params.workDir, 'uploads', fileName);
    if (!fs.existsSync(filePath)) continue;

    // Check by extension first
    let isFasta = isFastaFile(fileName);

    // If extension doesn't match, peek at content
    if (!isFasta) {
      try {
        const fd = fs.openSync(filePath, 'r');
        const buf = Buffer.alloc(512);
        const bytesRead = fs.readSync(fd, buf, 0, 512, 0);
        fs.closeSync(fd);
        const head = buf.toString('utf-8', 0, bytesRead);
        isFasta = looksLikeFasta(head);
      } catch {
        continue;
      }
    }

    if (!isFasta) continue;

    fastaFiles.push(fileName);

    // Count sequences (lines starting with ">")
    try {
      const content = fs.readFileSync(filePath, 'utf-8');
      const seqCount = (content.match(/^>/gm) ?? []).length;
      const lines = content.split('\n');

      // Detect if it's an alignment (all sequences same length after removing gaps info)
      const seqLengths = new Set<number>();
      let currentLen = 0;
      let inSeq = false;
      for (const line of lines) {
        if (line.startsWith('>')) {
          if (inSeq) seqLengths.add(currentLen);
          currentLen = 0;
          inSeq = true;
        } else {
          currentLen += line.trim().length;
        }
      }
      if (inSeq) seqLengths.add(currentLen);

      const isAligned = seqLengths.size === 1 && seqCount > 1;
      const ext = fileName.split('.').pop()?.toLowerCase() ?? '';
      const isAminoAcid = ext === 'faa' || fileName.includes('.faa');
      const stat = fs.statSync(filePath);
      const sizeKB = (stat.size / 1024).toFixed(1);

      const typeStr = isAligned ? 'aligned' : 'unaligned';
      const molStr = isAminoAcid ? 'amino acid' : 'nucleotide/unknown';

      detections.push(
        `- ${fileName}: FASTA (${typeStr}, ${molStr}), ${seqCount} sequences, ${sizeKB} KB` +
        (isAligned ? `, alignment length: ${[...seqLengths][0]} positions` : ''),
      );

      metadata[fileName] = {
        format: 'fasta',
        sequenceCount: seqCount,
        isAligned,
        isAminoAcid,
        alignmentLength: isAligned ? [...seqLengths][0] : null,
      };
    } catch {
      detections.push(`- ${fileName}: FASTA file (could not parse details)`);
    }
  }

  if (fastaFiles.length === 0) {
    return { contextText: '', metadata: {} };
  }

  return {
    contextText:
      `FASTA file detection:\n${detections.join('\n')}\n\n` +
      `NOTE: These are sequence files. Use Biopython (Bio.SeqIO / Bio.AlignIO) to parse them, NOT pandas. ` +
      `For alignment tasks, MAFFT is available as a CLI tool. For trimming, ClipKIT is installed as a Python package and CLI.`,
    metadata,
  };
}
