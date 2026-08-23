import { readFile, writeFile } from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-check-pattern-trace input.bin output.bin");
}
const source = await readFile(inputPath);
let offset = 0;
const take = (bytes) => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) {
    throw new Error(`truncated CheckPatternX trace at ${offset}`);
  }
  offset += bytes;
  return value;
};
const rawHeader = take(16);
if (rawHeader.readUInt32LE(0) !== 0x43504154 ||
    rawHeader.readUInt32LE(4) !== 1) {
  throw new Error("unexpected CheckPatternX trace header");
}
const nextNo = rawHeader.readInt32LE(8);
const sequenceLength = rawHeader.readInt32LE(12);
const count = nextNo + 1;
const positions = sequenceLength + 1;
const inputs = [
  take(12), take(12), take(12), take(2 * positions * count),
  take(8 * 9 * count), take(3 * count),
];
const outputs = [take(4), take(8 * 9 * count), take(3 * count)];
if (offset !== source.length) {
  throw new Error("unexpected trailing CheckPatternX trace data");
}
const header = Buffer.alloc(20);
header.write("CHKPAT1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeInt32LE(nextNo, 12);
header.writeInt32LE(sequenceLength, 16);
const sections = [];
const put = (id, value) => {
  const section = Buffer.alloc(8);
  section.writeUInt32LE(id, 0);
  section.writeUInt32LE(value.length, 4);
  sections.push(section, value);
};
inputs.forEach((value, index) => put(index + 1, value));
outputs.forEach((value, index) => put(index + 101, value));
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));
console.log(JSON.stringify({ nextNo, sequenceLength }));
