import { readFile, writeFile } from "node:fs/promises";

const [inputPath, outputPath, recordInputPath] = process.argv.slice(2);
if (!inputPath || !outputPath) {
  throw new Error("usage: convert-collectevents-boundary input.bin output.bin [records.bin]");
}
const source = await readFile(inputPath);
let offset = 0;
const take = (bytes) => {
  const value = source.subarray(offset, offset + bytes);
  if (value.length !== bytes) throw new Error(`truncated collectevents boundary at ${offset}`);
  offset += bytes;
  return value;
};
const records = [];
while (offset < source.length) {
  const rawHeader = take(17 * 4);
  if (rawHeader.readUInt32LE(0) !== 0x43454244) {
    throw new Error(`bad collectevents boundary magic at ${offset - 68}`);
  }
  const nextNo = rawHeader.readInt32LE(8);
  const sequenceLength = rawHeader.readInt32LE(12);
  const ubSmat = rawHeader.readInt32LE(32);
  const ubCollectEvents = rawHeader.readInt32LE(40);
  const count = nextNo + 1;
  const inputs = [
    take(6 * 4),
    take((sequenceLength + 1) * 4),
    take(6 * 4),
    take(3 * count * 4),
    take(3 * 4),
    take((ubSmat + 1) * count * 4),
    take(count * 2),
    take(count * 2),
    take((ubCollectEvents + 1) * 56),
  ];
  records.push({ rawHeader, nextNo, sequenceLength, inputs });
}
if (records.length < 2) throw new Error("expected two MakeCollecteventsC calls");
const selected = records.slice(0, 2);
if (selected.some((record) => record.nextNo !== selected[0].nextNo ||
    record.sequenceLength !== selected[0].sequenceLength)) {
  throw new Error("first MakeCollecteventsC calls have inconsistent dimensions");
}
const header = Buffer.alloc(20);
header.write("COLLECT1", 0, "ascii");
header.writeUInt32LE(1, 8);
header.writeInt32LE(selected[0].nextNo, 12);
header.writeInt32LE(selected[0].sequenceLength, 16);
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
  record.inputs.forEach((value, index) => putSection(base + index + 2, value));
}
if (recordInputPath) {
  const recordSource = await readFile(recordInputPath);
  if (recordSource.length < 52 || recordSource.readUInt32LE(0) !== 0x43434952) {
    throw new Error("bad collectevents input-record header");
  }
  const recordCount = recordSource.readUInt32LE(20);
  const recordBytes = recordSource.readUInt32LE(24);
  let recordOffset = 52;
  const locations = [];
  const eventRecords = [];
  const nonblankRdpCounts = new Int16Array(selected[0].nextNo + 1);
  for (let record = 0; record < recordCount; ++record) {
    if (recordOffset + 8 + recordBytes > recordSource.length) {
      throw new Error("truncated collectevents input record");
    }
    const location = recordSource.subarray(recordOffset, recordOffset + 8);
    const event = recordSource.subarray(
      recordOffset + 8, recordOffset + 8 + recordBytes);
    // Disabled-method calls in FinalTrim can increment CurrentXOver without
    // initializing the new PXOList slot.  The zero-filled slot subsequently
    // looks like ProgramFlag 0, but it is not an RDP event and cannot overlap
    // a real region.  Keep the RDP fixture semantic rather than counting that
    // unrelated native bookkeeping artifact as an event.
    const isBlank = event.every((value) => value === 0);
    if (event[2] === 0 && !isBlank) {
      locations.push(location);
      eventRecords.push(event);
      ++nonblankRdpCounts[location.readInt32LE(0)];
    }
    recordOffset += 8 + recordBytes;
  }
  putSection(11, Buffer.concat(locations));
  putSection(12, Buffer.concat(eventRecords));
  for (const record of selected) {
    for (let row = 0; row < nonblankRdpCounts.length; ++row) {
      record.inputs[7].writeInt16LE(nonblankRdpCounts[row], row * 2);
    }
  }
}
const end = Buffer.alloc(8);
end.writeUInt32LE(0xffffffff, 0);
sections.push(end);
await writeFile(outputPath, Buffer.concat([header, ...sections]));
console.log(JSON.stringify({
  nextNo: selected[0].nextNo,
  selected: [selected[0].rawHeader.readInt32LE(48),
    selected[0].rawHeader.readInt32LE(52),
    selected[0].rawHeader.readInt32LE(56)],
  parentRoles: selected.map((record) => record.rawHeader.readInt32LE(16)),
  rnum: [0, 1, 2].map((role) => selected[0].inputs[4].readInt32LE(role * 4)),
}));
