import { readFile, writeFile } from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-phpr-trace input.bin output.bin");
}

const source = await readFile(inputPath);
let offset = 0;
const take = (bytes) => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) {
    throw new Error(`truncated MakePhPrScore trace at ${offset}`);
  }
  offset += bytes;
  return value;
};

const records = [];
while (offset < source.length) {
  const rawHeader = take(16);
  if (rawHeader.readUInt32LE(0) !== 0x4d505052) {
    throw new Error(`bad MakePhPrScore magic at ${offset - 16}`);
  }
  const nextNo = rawHeader.readInt32LE(12);
  const count = nextNo + 1;
  records.push({
    rawHeader,
    nextNo,
    values: [
      take(8),
      take(3 * 4),
      take(2 * count * 4),
      take(count * 4),
      take(count * count * 4),
      take(count * count * 4),
      take(3 * 8),
      take(3 * 8),
      take(3 * 8),
    ],
  });
}
if (records.length < 2) {
  throw new Error(`expected at least two MakePhPrScore calls, got ${records.length}`);
}
const selected = records.slice(0, Math.min(3, records.length));
if (selected.some((record) => record.nextNo !== selected[0].nextNo)) {
  throw new Error("MakePhPrScore calls have inconsistent dimensions");
}

const header = Buffer.alloc(20);
header.write("PHPRSCO1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeInt32LE(selected[0].nextNo, 12);
header.writeUInt32LE(selected.length, 16);
const sections = [];
const putSection = (id, value) => {
  const section = Buffer.alloc(8);
  section.writeUInt32LE(id, 0);
  section.writeUInt32LE(value.length, 4);
  sections.push(section, value);
};
for (const [call, record] of selected.entries()) {
  const base = call * 1000;
  putSection(base + 1, record.rawHeader);
  record.values.forEach((value, index) => putSection(base + index + 2, value));
}
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));

console.log(JSON.stringify({
  nextNo: selected[0].nextNo,
  calls: selected.length,
  selected: selected[0].values[1].toString("hex").match(/.{8}/g)
    .map((hex) => Buffer.from(hex, "hex").readInt32LE()),
  scores: selected.map((record) => [0, 1, 2].map((index) =>
    record.values[6].readDoubleLE(index * 8))),
}));
