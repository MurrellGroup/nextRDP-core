import { readFile, writeFile } from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-modseqnumy-trace input.bin output.bin");
}
const source = await readFile(inputPath);
let offset = 0;
const take = (bytes) => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) throw new Error(`truncated ModSeqNumY trace at ${offset}`);
  offset += bytes;
  return value;
};
const calls = [];
while (offset < source.length) {
  const rawHeader = take(40);
  if (rawHeader.readUInt32LE(0) !== 0x4d535159) {
    throw new Error(`bad ModSeqNumY magic at ${offset - 40}`);
  }
  const length = rawHeader.readInt32LE(16);
  const role = rawHeader.readInt32LE(20);
  const last = [rawHeader.readInt32LE(24), rawHeader.readInt32LE(28),
    rawHeader.readInt32LE(32)];
  const inputCount = last.reduce((total, value) => total + value + 1, 0);
  const inputs = take(inputCount * 20);
  const result = take(4);
  const outputCount = last[role] + 1;
  const outputRecords = [];
  const sequenceRows = [];
  const savedRows = [];
  const missingRows = [];
  for (let record = 0; record < outputCount; ++record) {
    outputRecords.push(take(20));
    sequenceRows.push(take((length + 1) * 2));
    savedRows.push(take((length + 1) * 2));
    missingRows.push(take(length + 1));
  }
  calls.push({ rawHeader, inputs, result, outputRecords, sequenceRows,
    savedRows, missingRows });
}
const header = Buffer.alloc(20);
header.write("MSEQYV1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeUInt32LE(calls.length, 12);
const sections = [];
const putSection = (id, value) => {
  const section = Buffer.alloc(8);
  section.writeUInt32LE(id, 0);
  section.writeUInt32LE(value.length, 4);
  sections.push(section, value);
};
for (const [call, record] of calls.entries()) {
  const base = call * 1000;
  putSection(base + 1, record.rawHeader);
  putSection(base + 2, record.inputs);
  putSection(base + 3, record.result);
  putSection(base + 4, Buffer.concat(record.outputRecords));
  putSection(base + 5, Buffer.concat(record.sequenceRows));
  putSection(base + 6, Buffer.concat(record.savedRows));
  putSection(base + 7, Buffer.concat(record.missingRows));
}
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));
console.log(JSON.stringify({ calls: calls.length }));
