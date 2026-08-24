import { readFile, writeFile } from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-mutation-tail-trace input.bin output.bin");
}
const source = await readFile(inputPath);
let offset = 0;
const take = (bytes) => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) throw new Error(`truncated mutation tail at ${offset}`);
  offset += bytes;
  return value;
};
const events = [];
while (offset < source.length) {
  const msnHeader = take(44);
  if (msnHeader.readUInt32LE(0) !== 0x4d534e5a) {
    throw new Error(`expected ModSN at ${offset - 44}`);
  }
  const length = msnHeader.readInt32LE(12);
  const role = msnHeader.readInt32LE(24);
  const count = msnHeader.readInt32LE(28 + role * 4) + 1;
  const msnRecords = [], newSequences = [], newMissing = [];
  for (let record = 0; record < count; ++record) {
    msnRecords.push(take(20));
    newSequences.push(take((length + 1) * 2));
    newMissing.push(take(length + 1));
  }
  const zHeader = take(44);
  if (zHeader.readUInt32LE(0) !== 0x4d53515a) {
    throw new Error(`expected ModSeqNumZ at ${offset - 44}`);
  }
  const zRecords = [], originalSequences = [], originalMissing = [];
  for (let record = 0; record < count; ++record) {
    zRecords.push(take(16));
    originalSequences.push(take((length + 1) * 2));
    originalMissing.push(take(length + 1));
  }
  const sizeHeader = take(36);
  if (sizeHeader.readUInt32LE(0) !== 0x4d415353) {
    throw new Error(`expected MakeActualSeqSize at ${offset - 36}`);
  }
  const sizeRecords = take(count * 20);
  events.push({ msnHeader, msnRecords, newSequences, newMissing, zHeader,
    zRecords, originalSequences, originalMissing, sizeHeader, sizeRecords });
}
const header = Buffer.alloc(20);
header.write("MTAILV1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeUInt32LE(events.length, 12);
const sections = [];
const put = (id, value) => {
  const section = Buffer.alloc(8);
  section.writeUInt32LE(id, 0);
  section.writeUInt32LE(value.length, 4);
  sections.push(section, value);
};
for (const [event, record] of events.entries()) {
  const base = event * 1000;
  put(base + 1, record.msnHeader);
  put(base + 2, Buffer.concat(record.msnRecords));
  put(base + 3, Buffer.concat(record.newSequences));
  put(base + 4, Buffer.concat(record.newMissing));
  put(base + 5, record.zHeader);
  put(base + 6, Buffer.concat(record.zRecords));
  put(base + 7, Buffer.concat(record.originalSequences));
  put(base + 8, Buffer.concat(record.originalMissing));
  put(base + 9, record.sizeHeader);
  put(base + 10, record.sizeRecords);
}
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));
console.log(JSON.stringify({ events: events.length }));
