#!/usr/bin/env node

import fs from "node:fs";

const path = process.argv[2];
if (!path) {
  console.error("usage: inspect-make-test-pvs.mjs FIXTURE");
  process.exit(2);
}

const data = fs.readFileSync(path);
if (data.subarray(0, 8).toString("ascii") !== "MKTESTPV") {
  throw new Error("not a MakeTestPVs fixture");
}
const version = data.readUInt32LE(8);
const doneUb = data.readInt32LE(12);
const nextNo = data.readInt32LE(16);
const rowUb = data.readInt32LE(20);
const slotUb = data.readInt32LE(24);
const eventBytes = data.readInt32LE(28);
if (version !== 1 || eventBytes !== 56) {
  throw new Error(`unsupported fixture version/record size ${version}/${eventBytes}`);
}

const sections = new Map();
let offset = 32;
while (offset + 8 <= data.length) {
  const id = data.readUInt32LE(offset);
  const bytes = data.readUInt32LE(offset + 4);
  offset += 8;
  if (id === 0xffffffff) break;
  if (offset + bytes > data.length) throw new Error("truncated section");
  sections.set(id, data.subarray(offset, offset + bytes));
  offset += bytes;
}

const current = sections.get(2);
const events = sections.get(3);
if (!current || !events) throw new Error("fixture has no input event state");
const rowCount = rowUb + 1;
const totals = new Map();
let total = 0;
console.log(`nextNo=${nextNo} rowUb=${rowUb} slotUb=${slotUb} doneUb=${doneUb}`);
for (let row = 0; row <= nextNo; ++row) {
  const count = current.readInt16LE(row * 2);
  total += count;
  if (count) console.log(`row ${row}: ${count}`);
  for (let slot = 1; slot <= count; ++slot) {
    const base = (row + slot * rowCount) * eventBytes;
    const program = events.readUInt8(base + 2);
    totals.set(program, (totals.get(program) ?? 0) + 1);
    const major = events.readInt16LE(base + 6);
    const minor = events.readInt16LE(base + 8);
    const daughter = events.readInt16LE(base + 10);
    const beginning = events.readInt32LE(base + 12);
    const ending = events.readInt32LE(base + 16);
    const probability = events.readDoubleLE(base + 40);
    console.log(
      `  ${slot}: p${program} d=${daughter} mi=${minor} ma=${major} ` +
      `b=${beginning} e=${ending} prob=${probability}`,
    );
  }
}
console.log(`total=${total}`);
console.log(
  `programs=${[...totals.entries()].sort((a, b) => a[0] - b[0])
    .map(([program, count]) => `${program}:${count}`).join(",")}`,
);
