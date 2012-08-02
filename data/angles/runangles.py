#!/usr/bin/regina-python --nolibs

import sys

if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0] + " <filename>"
    sys.exit(1)

tree = regina.readFileMagic(sys.argv[1])
if not tree:
    print "E: Could not open file " + sys.argv[1] + "."
    print
    print "Usage: " + sys.argv[0] + " <filename>"
    sys.exit(1)

def process(tri):
    a = regina.NAngleStructureList.enumerate(tri)
    list = []
    for i in range(a.getNumberOfStructures()):
        list.append(a.getStructure(i).toStringLong())
    list.sort()
    for i in list:
        sys.stdout.write(i)
    sys.stdout.flush()

p = tree
while p != None:
    if p.getPacketType() == regina.NTriangulation.packetType:
        print p.getPacketLabel()

        print "Angle structures:"
        process(p)

    p = p.nextTreePacket()

