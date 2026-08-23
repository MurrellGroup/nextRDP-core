import {readFile, writeFile} from "node:fs/promises";

const [inputPath, outputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-consensus-csv NN_inputs.txt output.bin");
}
const lines = (await readFile(inputPath, "utf8"))
  .split(/\r?\n/)
  .filter(line => line.trim().length > 0);
if (lines.length < 2) throw new Error("consensus CSV has no data row");
const values = lines[1].split(",").map(value => Number(value.trim()));
if (values.length !== 120 || values.some(value => !Number.isFinite(value))) {
  throw new Error(`expected 120 finite first-event values, found ${values.length}`);
}
const header = Buffer.alloc(12);
header.write("CONSCV1", 0, "ascii");
header.writeUInt32LE(1, 8);
const payload = Buffer.alloc(values.length * 8);
values.forEach((value, index) => payload.writeDoubleLE(value, index * 8));
const section = Buffer.alloc(8);
section.writeUInt32LE(1, 0);
section.writeUInt32LE(payload.length, 4);
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
await writeFile(outputPath, Buffer.concat([header, section, payload, end]));
console.log(JSON.stringify({values: values.length}));
