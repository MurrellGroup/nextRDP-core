import { readFile, writeFile } from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-strip-dup-inv-trace input.bin output.bin");
}
const source = await readFile(inputPath);
if (source.length < 16 || source.readUInt32LE(0) !== 0x53444956) {
  throw new Error("input does not begin with an RDP StripDupInv trace record");
}
const nextNo = source.readInt32LE(12);
const count = nextNo + 1;
let offset = 16;
const take = (bytes) => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) throw new Error("truncated StripDupInv trace");
  offset += bytes;
  return value;
};
const inputs = [take(12), take(12 * count), take(12 * count)];
const outputs = [take(12), take(12 * count), take(12 * count), take(12)];

const header = Buffer.alloc(16);
header.write("STRIPDI1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeInt32LE(nextNo, 12);
const sections = [];
const putSection = (id, value) => {
  const section = Buffer.alloc(8);
  section.writeUInt32LE(id, 0);
  section.writeUInt32LE(value.length, 4);
  sections.push(section, value);
};
inputs.forEach((value, index) => putSection(index + 1, value));
outputs.forEach((value, index) => putSection(index + 101, value));
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));
