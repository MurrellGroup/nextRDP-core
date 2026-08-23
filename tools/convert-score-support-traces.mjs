import { readFile, writeFile } from "node:fs/promises";

const [donePath, groupPath, scorePath, outputPath] = process.argv.slice(2);
if (!donePath || !groupPath || !scorePath || !outputPath) {
  throw new Error("usage: convert-score-support-traces done.bin groups.bin score.bin output.bin");
}

const parse = async (path, magic, readRecord) => {
  const source = await readFile(path);
  let offset = 0;
  const take = (bytes) => {
    const value = source.subarray(offset, offset + bytes);
    if (value.length !== bytes) throw new Error(`truncated trace ${path} at ${offset}`);
    offset += bytes;
    return value;
  };
  const records = [];
  while (offset < source.length) {
    const rawHeader = take(16);
    if (rawHeader.readUInt32LE(0) !== magic) {
      throw new Error(`bad trace magic in ${path} at ${offset - 16}`);
    }
    records.push({rawHeader, ...readRecord(rawHeader, take)});
  }
  return records;
};

const doneRecords = await parse(donePath, 0x4d445433, (header, take) => {
  const nextNo = header.readInt32LE(8);
  const count = nextNo + 1;
  return {nextNo, values: [take(12), take(2 * count * 4),
    take(3 * count * 4), take(3 * count * 4),
    take(3 * count * 4), take(4), take(2 * count * 4)]};
});
const groupRecords = await parse(groupPath, 0x54524732, (header, take) => {
  const nextNo = header.readInt32LE(12);
  const count = nextNo + 1;
  return {nextNo, values: [take(12), take(24), take(count * 4),
    take(count * 4), take(count * 4), take(24),
    take(3 * count * 4), take(4), take(count * 4),
    take(count * 4), take(count * 4), take(24)]};
});
const scoreRecords = await parse(scorePath, 0x54525332, (header, take) => {
  const nextNo = header.readInt32LE(12);
  const count = nextNo + 1;
  return {nextNo, values: [take(12), take(32), take(count * 4),
    take(count * 4), take(3 * count * 4), take(3 * count * 4),
    take(4), take(32)]};
});

const done = doneRecords.slice(0, 2);
const groups = groupRecords.slice(0, 3);
const scores = scoreRecords.slice(0, 3);
if (done.length !== 2 || groups.length !== 3 || scores.length !== 3) {
  throw new Error(`unexpected score support calls: ${done.length}/${groups.length}/${scores.length}`);
}
const nextNo = done[0].nextNo;
if ([...done, ...groups, ...scores].some((record) => record.nextNo !== nextNo)) {
  throw new Error("score support calls have inconsistent dimensions");
}

const header = Buffer.alloc(28);
header.write("SCORESP1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeInt32LE(nextNo, 12);
header.writeUInt32LE(done.length, 16);
header.writeUInt32LE(groups.length, 20);
header.writeUInt32LE(scores.length, 24);
const sections = [];
const putSection = (id, value) => {
  const section = Buffer.alloc(8);
  section.writeUInt32LE(id, 0);
  section.writeUInt32LE(value.length, 4);
  sections.push(section, value);
};
const putRecords = (records, familyBase) => records.forEach((record, call) => {
  const base = familyBase + call * 1000;
  putSection(base + 1, record.rawHeader);
  record.values.forEach((value, index) => putSection(base + index + 2, value));
});
putRecords(done, 0);
putRecords(groups, 10000);
putRecords(scores, 20000);
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));
console.log(JSON.stringify({nextNo, done: done.length, groups: groups.length,
  scores: scores.length, selected: [0, 4, 8].map((offset) =>
    done[0].values[0].readInt32LE(offset))}));
