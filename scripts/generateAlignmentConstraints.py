import sys
import io
import numpy as np

def get_alignconsts(filename, tablename):
    iobuf = io.StringIO()
    read_lines = False
    with open(filename, 'r') as f:
        for line in f.readlines():
            if 'TABLE' in line:
                if line.strip().split()[-1] == tablename:
                    read_lines = True
                    continue
                else:
                    read_lines = False
                    continue

            if read_lines:
              iobuf.write(','.join(line.split(",")[2:]))
    iobuf.seek(0)

    return np.genfromtxt(iobuf,delimiter=',')

fn = sys.argv[1]
plane_alignment = get_alignconsts(fn,"TrkAlignPlane")
panel_alignment = get_alignconsts(fn,"TrkAlignPanel")

xskew = plane_alignment[0][0]-plane_alignment[-1][0]
yskew = plane_alignment[0][1]-plane_alignment[-1][1]
zsqueeze = plane_alignment[0][2]-plane_alignment[-1][2]
ztwist = plane_alignment[0][5]-plane_alignment[-1][5]

rsqueeze = np.sum(panel_alignment[:,1])

fout = open("constraints.fcl","w")
fout.write("physics.analyzers.AlignTrackCollector.WeakValues : [\n")
fout.write("  %f, %f, %f, %f, %f]\n" % (xskew, yskew, zsqueeze, ztwist, rsqueeze))
fout.write("physics.analyzers.AlignTrackCollector.PanelValuesX : [\n")
for i in range(len(panel_alignment[:,0])):
  fout.write(" %f" % panel_alignment[i][0])
  if i != len(panel_alignment[:,0])-1:
    fout.write(",")
fout.write("]\n")

fout.write("physics.analyzers.AlignTrackCollector.PanelValuesY : [\n")
for i in range(len(panel_alignment[:,0])):
  fout.write(" %f" % panel_alignment[i][1])
  if i != len(panel_alignment[:,0])-1:
    fout.write(",")
fout.write("]\n")

fout.write("physics.analyzers.AlignTrackCollector.PanelValuesZ : [\n")
for i in range(len(panel_alignment[:,0])):
  fout.write(" %f" % panel_alignment[i][2])
  if i != len(panel_alignment[:,0])-1:
    fout.write(",")
fout.write("]\n")

fout.write("physics.analyzers.AlignTrackCollector.PanelValuesRX : [\n")
for i in range(len(panel_alignment[:,0])):
  fout.write(" %f" % panel_alignment[i][3])
  if i != len(panel_alignment[:,0])-1:
    fout.write(",")
fout.write("]\n")

fout.write("physics.analyzers.AlignTrackCollector.PanelValuesRY : [\n")
for i in range(len(panel_alignment[:,0])):
  fout.write(" %f" % panel_alignment[i][4])
  if i != len(panel_alignment[:,0])-1:
    fout.write(",")
fout.write("]\n")

fout.write("physics.analyzers.AlignTrackCollector.PanelValuesRZ : [\n")
for i in range(len(panel_alignment[:,0])):
  fout.write(" %f" % panel_alignment[i][5])
  if i != len(panel_alignment[:,0])-1:
    fout.write(",")
fout.write("]\n")

fout.close()
