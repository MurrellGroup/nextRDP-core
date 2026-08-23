import {readFile, writeFile} from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-cmaxd2p3-trace input.bin output.bin");
}
const source = await readFile(inputPath);
let offset = 0;
const take = bytes => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) throw new Error(`truncated CMaxD2P3 trace at ${offset}`);
  offset += bytes;
  return value;
};
const rawHeader = take(10 * 4);
if (rawHeader.readUInt32LE(0) !== 0x434d4433) {
  throw new Error("bad CMaxD2P3 trace magic");
}
const incnum = rawHeader.readInt32LE(8);
const nextNo = rawHeader.readInt32LE(32);
const sequenceLength = rawHeader.readInt32LE(36);
const count = nextNo + 1;
const sections = [
  take(count * (sequenceLength + 1) * 2),
  take((sequenceLength + 1) * 4),
  take((sequenceLength + 1) * 4),
  take(86),
  take((incnum + 1) * 4),
  take(count),
  take(count),
  take(1875 * 4),
  take(4),
  take(3 * 4),
  take(3 * 4),
];
const header = Buffer.alloc(20);
header.write("CMAXD3V1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeInt32LE(nextNo, 12);
header.writeInt32LE(sequenceLength, 16);
const output = [header];
const put = (id, value) => {
  const sectionHeader = Buffer.alloc(8);
  sectionHeader.writeUInt32LE(id, 0);
  sectionHeader.writeUInt32LE(value.length, 4);
  output.push(sectionHeader, value);
};
put(1, rawHeader);
sections.forEach((value, index) => put(index + 2, value));
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
output.push(end);
await writeFile(outputPath, Buffer.concat(output));
console.log(JSON.stringify({nextNo, sequenceLength, incnum, bytes: offset}));
