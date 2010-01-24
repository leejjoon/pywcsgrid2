import glob
import os.path

src_dir = "axes_grid_orig"
dest_dir = "axes_grid"

for fn in glob.glob(os.path.join(src_dir, "*.py")):
    #print fn
    basename = os.path.basename(fn)
    s = open(os.path.join(src_dir, basename)).read()
    s2 = s.replace("mpl_toolkits.axes_grid", "pywcsgrid2.axes_grid")
    open(os.path.join(dest_dir, basename),"w").write(s2)
