import { readFile, writeFile } from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-rcompat-trace input.bin output.bin");
}
const source = await readFile(inputPath);
let offset = 0;
const take = (bytes) => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) throw new Error(`truncated MakeRCompat trace at ${offset}`);
  offset += bytes;
  return value;
};
const records = [];
while (offset < source.length) {
  const rawHeader = take(16);
  if (rawHeader.readUInt32LE(0) !== 0x52434d50) {
    throw new Error(`bad MakeRCompat magic at ${offset - 16}`);
  }
  const invocation = rawHeader.readUInt32LE(4);
  const role = rawHeader.readInt32LE(8);
  const nextNo = rawHeader.readInt32LE(12);
  const count = nextNo + 1;
  const inputs = [
    take(12), take(24), take(12), take(12), take(12), take(12), take(12),
    take(8 * count), take(4 * count), take(12 * count), take(12 * count),
    take(4 * count * count), take(24),
  ];
  const outputs = [take(12), take(12), take(12)];
  records.push({ invocation, role, nextNo, inputs, outputs });
}
const firstSequences = records[0]?.inputs[0];
const boundary = records.findIndex((record, index) => index > 0 &&
  (record.nextNo !== records[0].nextNo ||
   !record.inputs[0].equals(firstSequences)));
const selectedRecords = records.slice(0, boundary < 0 ? records.length : boundary);
if (selectedRecords.length < 6 || selectedRecords.some((record, index) =>
    record.invocation !== index + 1 ||
    record.nextNo !== selectedRecords[0].nextNo)) {
  throw new Error("unexpected MakeRCompat call sequence");
}
const header = Buffer.alloc(20);
header.write("RCOMPAT1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeInt32LE(selectedRecords[0].nextNo, 12);
header.writeUInt32LE(selectedRecords.length, 16);
const sections = [];
const putSection = (id, value) => {
  const section = Buffer.alloc(8);
  section.writeUInt32LE(id, 0);
  section.writeUInt32LE(value.length, 4);
  sections.push(section, value);
};
for (const [call, record] of selectedRecords.entries()) {
  const base = call * 1000;
  const metadata = Buffer.alloc(8);
  metadata.writeInt32LE(record.invocation, 0);
  metadata.writeInt32LE(record.role, 4);
  putSection(base + 1, metadata);
  record.inputs.forEach((value, index) => putSection(base + index + 2, value));
  record.outputs.forEach((value, index) => putSection(base + index + 101, value));
}
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));
console.log(JSON.stringify({ nextNo: selectedRecords[0].nextNo,
  calls: selectedRecords.length,
  roles: selectedRecords.map((record) => record.role) }));
