import fs from 'node:fs';
import path from 'node:path';

const dir = path.dirname(new URL(import.meta.url).pathname.replace(/^\//, ''));
const read = (f) => fs.readFileSync(path.join(dir, f), 'utf8');

function parseFasta(text) {
  const out = [];
  let cur = null;
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line) continue;
    if (line.startsWith('>')) {
      const body = line.slice(1).trim();
      const sp = body.indexOf(' ');
      cur = { name: sp > 0 ? body.slice(0, sp) : body, desc: sp > 0 ? body.slice(sp + 1) : '', seq: '' };
      out.push(cur);
    } else if (cur) {
      cur.seq += line.replace(/[^A-Za-z]/g, '');
    }
  }
  return out;
}

const puc = parseFasta(read('L09137.fasta'))[0];
const phi = parseFasta(read('phiX174.fasta'))[0];
const ins = parseFasta(read('ins.fasta'))[0];

const mcs = 'GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT';
const norm = (s) => s.toUpperCase().replace(/[^ACGTNRYKMSWBDHV]/g, '');

const samples = [
  {
    id: 'puc19',
    name: 'pUC19',
    label: 'pUC19 克隆载体',
    desc: 'Cloning vector pUC19, complete sequence (ENA L09137.2)',
    topology: 'circular',
    seq: norm(puc.seq),
  },
  {
    id: 'phix174',
    name: 'phiX174',
    label: 'φX174 噬菌体基因组',
    desc: 'Coliphage phi-X174, complete genome (ENA J02482.1)',
    topology: 'circular',
    seq: norm(phi.seq),
  },
  {
    id: 'ins',
    name: 'INS',
    label: '人 胰岛素 mRNA (INS)',
    desc: 'Homo sapiens insulin, transcript ENST00000381330.5 (cDNA)',
    topology: 'linear',
    seq: norm(ins.seq),
  },
  {
    id: 'mcs',
    name: 'MCS',
    label: '多克隆位点片段 (MCS)',
    desc: 'pUC 系列多克隆位点 EcoRI–HindIII 片段（10 个常见酶切位点）',
    topology: 'linear',
    seq: mcs,
  },
];

for (const s of samples) console.log(`${s.id}: ${s.seq.length} bp, ${s.topology}`);

const body = samples
  .map((s) => {
    const chunks = s.seq.match(/.{1,100}/g).map((c) => `'${c}'`).join(' +\n    ');
    return `  {
    id: '${s.id}',
    name: '${s.name}',
    label: '${s.label}',
    desc: '${s.desc}',
    topology: '${s.topology}',
    seq:
    ${chunks},
  },`;
  })
  .join('\n');

const out = `/* 示例序列：由 _data/gen-samples.mjs 生成，请勿手工编辑。
   全部来自公开数据库（ENA pUC19=L09137.2 / phiX174=J02482.1 / Ensembl INS=ENST00000381330.5）。 */
const SAMPLE_SEQUENCES = [
${body}
];
`;

const target = path.resolve(dir, '..', '01-SeqLab', 'samples.js');
fs.writeFileSync(target, out, 'utf8');
console.log('written:', target);
