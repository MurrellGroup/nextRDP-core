import { readFile, writeFile } from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-find-actual-events-trace input.bin output.bin");
}

const source = await readFile(inputPath);
let offset = 0;
const records = [];
const take = (bytes) => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) {
    throw new Error(`truncated FindActualEvents trace at ${offset}`);
  }
  offset += bytes;
  return value;
};

while (offset < source.length) {
  const header = take(24);
  if (header.readUInt32LE(0) !== 0x46414556) {
    throw new Error(`bad FindActualEvents trace magic at ${offset - 24}`);
  }
  const invocation = header.readUInt32LE(4);
  const sequenceLength = header.readInt32LE(8);
  const winPp = header.readInt32LE(12);
  const nextNo = header.readInt32LE(16);
  const ub = header.readInt32LE(20);
  const count = nextNo + 1;
  const positions = sequenceLength + 1;
  const inputs = [
    take(3 * 4),
    take(6 * 4),
    take(6 * 4),
    take(3 * count),
    take(6 * count * 4),
    take(3 * count * 4),
    take(57 * count * 8),
    take(count * 4),
    take(6 * 4),
    take(6 * 4),
    take(9 * count * 4),
    take(positions * 4),
    take(positions * 4),
    take(positions * 4),
    take(2 * 4),
    take(3 * 4),
    take(3 * count * 4),
    take(3 * count),
    take(3 * 4),
    take(2 * 8),
    take(count * 2),
    take(4 * 4),
    take(6),
    take(3 * count),
  ];
  const eventCount = take(4).readUInt32LE(0);
  const events = take(eventCount * 64);
  const outputs = [
    take(4),
    take(6 * 4),
    take(6 * count * 4),
    take(3 * count * 4),
    take(57 * count * 8),
    take(count * 4),
    take(2 * 4),
    take(3 * 4),
    take(2 * 8),
    take(4 * 4),
    take(6),
  ];
  records.push({
    invocation,
    sequenceLength,
    winPp,
    nextNo,
    ub,
    eventCount,
    inputs,
    events,
    outputs,
  });
}

if (records.length !== 3 || records.some((record, index) =>
    record.invocation !== index + 1 || record.winPp !== index)) {
  throw new Error("expected the first three role calls in the trace");
}
const first = records[0];
if (records.some((record) =>
    record.sequenceLength !== first.sequenceLength ||
    record.nextNo !== first.nextNo || record.ub !== first.ub)) {
  throw new Error("FindActualEvents trace calls have inconsistent dimensions");
}

const header = Buffer.alloc(24);
header.write("FAEVENT1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeInt32LE(first.sequenceLength, 12);
header.writeInt32LE(first.nextNo, 16);
header.writeInt32LE(first.ub, 20);
const sections = [];
const putSection = (id, value) => {
  const section = Buffer.alloc(8);
  section.writeUInt32LE(id, 0);
  section.writeUInt32LE(value.length, 4);
  sections.push(section, value);
};
for (const [call, record] of records.entries()) {
  const base = call * 1000;
  record.inputs.forEach((value, index) => putSection(base + index + 1, value));
  const eventHeader = Buffer.alloc(4);
  eventHeader.writeUInt32LE(record.eventCount, 0);
  putSection(base + 25, Buffer.concat([eventHeader, record.events]));
  record.outputs.forEach((value, index) =>
    putSection(base + 101 + index, value));
}
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));

const readInts = (value) => Array.from(
  { length: value.length / 4 }, (_, index) => value.readInt32LE(index * 4));
for (const record of records) {
  const found = readInts(record.outputs[5]);
  console.log(JSON.stringify({
    invocation: record.invocation,
    role: record.winPp,
    selected: readInts(record.inputs[0]),
    rnum: readInts(record.inputs[15]),
    eventCount: record.eventCount,
    foundSlots: found.flatMap((value, index) => value ? [index] : []),
  }));
}
