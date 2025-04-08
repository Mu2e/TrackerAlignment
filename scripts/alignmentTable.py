import numpy as np
import sys

class AlignmentSet:
    def __init__(self):
        self.tables = {}
        self.xskew = 0
        self.yskew = 0
        self.xparallel = 0
        self.yparallel = 0
        self.zsqueeze = 0
        self.twist = 0
        self.xshiftrms = 0
        self.yshiftrms = 0
        self.zshiftrms = 0
        self.xanglerms = 0
        self.yanglerms = 0
        self.zanglerms = 0
        self.vshiftrms = 0
        self.wshiftrms = 0
        self.uanglerms = 0
        self.vanglerms = 0
        self.wanglerms = 0
        self.wirevrms = 0
        self.wirewrms = 0
        self.strawvrms = 0
        self.strawwrms = 0

    def fromFile(self,filename): 
        tablenames = ["TrkAlignTracker","TrkAlignPlane","TrkAlignPanel","TrkAlignStraw"]
        constants = parse_tables(filename,tablenames)
        tables = {}
        for name in tablenames:
            if name == "TrkAlignTracker":
                class_id = 0
                objects = 1
                ndof = 6
                stride = 0
            elif name == "TrkAlignPlane":
                class_id = 1
                objects = 36
                ndof = 6
                stride = 6*96
            elif name == "TrkAlignPanel":
                class_id = 2
                objects = 216
                ndof = 6
                stride = 96
            elif name == "TrkAlignStraw":
                class_id = 3
                objects = 20736
                ndof = 8
                stride = 1
            else:
                continue
            if name in constants:
                self.tables[name] = AlignmentTable(name,class_id,objects,ndof,stride,constants[name])

    def getStats(self):
        if "TrkAlignPlane" in self.tables:
            plane_data = self.tables["TrkAlignPlane"].constants.copy()
            for i in [0,2,3,5]:
                plane_data[:,i][1::4] *= -1
                plane_data[:,i][2::4] *= -1

            lm = np.polyfit(np.arange(0,len(plane_data[:,0]),1),plane_data[:,0],1)
            fn = np.poly1d(lm)
            self.xskew = fn(len(plane_data[:,0])-1)-fn(0)
            plane_data[:,0] -= fn(np.arange(0,len(plane_data[:,0]),1))
            lm = np.polyfit(np.arange(0,len(plane_data[:,1]),1),plane_data[:,1],1)
            fn = np.poly1d(lm)
            self.yskew = fn(len(plane_data[:,1])-1)-fn(0)
            plane_data[:,1] -= fn(np.arange(0,len(plane_data[:,1]),1))
            lm = np.polyfit(np.arange(0,len(plane_data[:,2]),1),plane_data[:,2],1)
            fn = np.poly1d(lm)
            self.zsqueeze = fn(len(plane_data[:,2])-1)-fn(0)
            plane_data[:,2] -= fn(np.arange(0,len(plane_data[:,2]),1))
            lm = np.polyfit(np.arange(0,len(plane_data[:,5]),1),plane_data[:,5],1)
            fn = np.poly1d(lm)
            plane_data[:,5] -= fn(np.arange(0,len(plane_data[:,5]),1))
            self.twist = fn(len(plane_data[:,5])-1)-fn(0)
            self.xparallel = np.mean(plane_data[:,3])
            self.yparallel = np.mean(plane_data[:,4])
            
            self.xshiftrms = np.std(plane_data[:,0])
            self.yshiftrms = np.std(plane_data[:,1])
            self.zshiftrms = np.std(plane_data[:,2])
            self.xanglerms = np.std(plane_data[:,3])
            self.yanglerms = np.std(plane_data[:,4])
            self.zanglerms = np.std(plane_data[:,5])

        if "TrkAlignPanel" in self.tables:
            panel_data = self.tables["TrkAlignPanel"].constants.copy()
            self.rsqueeze = np.mean(panel_data[:,1])
            self.vshiftrms = np.std(panel_data[:,1])
            self.wshiftrms = np.std(panel_data[:,2])
            self.uanglerms = np.std(panel_data[:,3])
            self.vanglerms = np.std(panel_data[:,4])
            self.wanglerms = np.std(panel_data[:,5])

        if "TrkAlignStraw" in self.tables:
            self.strawvrms = np.std(np.concatenate([straw_data[:,4],straw_data[:,6]]))
            self.strawvrms = np.std(np.concatenate([straw_data[:,5],straw_data[:,7]]))
            if "TrkAlignStrawSim" in self.tables:
                # FIXME
                self.wirevrms = np.std(np.concatenate([straw_data[:,0],straw_data[:,2]]))
                self.wirewrms = np.std(np.concatenate([straw_data[:,1],straw_data[:,3]]))

    def printStats(self):
        self.getStats()
        if "TrkAlignPlane" in self.tables:
            print("PLANE")
            print("  x translation rms (um): %7.1f" % (self.xshiftrms*1e3))
            print("  y translation rms (um): %7.1f" % (self.yshiftrms*1e3))
            print("  z translation rms (um): %7.1f" % (self.zshiftrms*1e3))
            print("  x rotation rms    (ur): %7.1f" % (self.xanglerms*1e6))
            print("  y rotation rms    (ur): %7.1f" % (self.yanglerms*1e6))
            print("  z rotation rms    (ur): %7.1f" % (self.zanglerms*1e6))
            print("PLANE WEAK MODES")
            print("  x skew            (um): %7.1f" % (self.xskew*1e3))
            print("  y skew            (um): %7.1f" % (self.yskew*1e3))
            print("  z squeeze         (um): %7.1f" % (self.zsqueeze*1e3))
            print("  x parallel        (ur): %7.1f" % (self.xparallel*1e6))
            print("  y parallel        (ur): %7.1f" % (self.yparallel*1e6))
            print("  z twist           (ur): %7.1f" % (self.twist*1e6))
        if "TrkAlignPanel" in self.tables:
            print("PANEL")
            print("  v translation rms (um): %7.1f" % (self.vshiftrms*1e3))
            print("  w translation rms (um): %7.1f" % (self.wshiftrms*1e3))
            print("  u rotation rms    (ur): %7.1f" % (self.uanglerms*1e6))
            print("  v rotation rms    (ur): %7.1f" % (self.vanglerms*1e6))
            print("  w rotation rms    (ur): %7.1f" % (self.wanglerms*1e6))
            print("PANEL WEAK MODES")
            print("  r sqeeze          (um): %7.1f" % (self.rsqueeze*1e3))
        if "TrkAlignStraw" in self.tables:
            print("STRAW")
            print("  straw v rms       (um): %7.1f" % (self.strawvrms*1e3))
            print("  straw w rms       (um): %7.1f" % (self.strawwrms*1e3))
            if "TrkAlignStrawSim" in self.tables:
                print("  wire v rms        (um): %7.1f" % (self.wirevrms*1e3))
                print("  wire w rms        (um): %7.1f" % (self.wirewrms*1e3))

class AlignmentTable:
    def __init__(self, table, class_id, objects, ndof, stride, constants=None, errors=None):
        self.table = table
        self.nobjects = objects
        self.ndof = ndof
        self.classid = class_id
        self.stride = stride

        self.constants = []
        self.errors = []
        if not constants is None:
            self.constants = constants
        else:
            for _ in range(self.nobjects):
                self.constants.append([0.0]*ndof)
        if not errors is None:
            self.errors = errors
        else:
            for _ in range(self.nobjects):
                self.errors.append([0.0]*ndof)

        self.constants = np.array(self.constants)
        self.errors = np.array(self.errors)

    def table_name(self):
        return self.table 
    
    def n_objects(self):
        return self.nobjects
    
    def dof_per_obj(self):
        return self.ndof
    
    def mplabel(self, id, dof):
        return 10000 * self.classid + 10 * id + dof

    def strawId(self, index):
        uniquestraw = index*self.stride
        plane = int(uniquestraw/(6*96))
        panel = int(uniquestraw/96)%6
        straw = uniquestraw%96
        return "%d_%d_%d" % (plane,panel,straw)
    
    def setConstant(self, id, dof, value):
        self.constants[id][dof] = float(value)
    
    def setError(self, id, dof, error):
        self.errors[id][dof] = float(error)

    def to_proditions_table(self):
        lines = []
        for i, constants in enumerate(self.constants):
            lines.append('%d,' % i + self.strawId(i) + ',' +  ','.join([str(i) for i in constants]))
        return """TABLE {name}
{csv}

""".format(name=self.table,
            csv='\n'.join(lines))

def parse_tables(filename, table_names=None):
    result = {}
    current_table = None
    rows = []
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            
            if not line or line.startswith('#'):
                continue
                
            if line.startswith('TABLE '):
                if current_table and rows:
                    if current_table in table_names:
                        result[current_table] = np.array(rows,dtype=float)
                
                current_table = line.strip().split()[1]
                rows = []
            else:
                parts = line.split(',')
                data_row = parts[2:]
                rows.append(data_row)
    
    if current_table and rows:
        if current_table in table_names:
            result[current_table] = np.array(rows,dtype=float)
    
    return result

if __name__ == "__main__":
    aligns = AlignmentSet()
    for fn in sys.argv[1:]:
        aligns.fromFile(fn)
    aligns.printStats()
