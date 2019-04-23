import glob, os, subprocess

# tag = "test_v2"
# stlTag = "TDR-3x2_1190mmOuterRad"
tag = "barrelTracks_r1p16_eta1p8_pt0p7"
stlTag = "BTL-1Tray_1160mmOuterRad"

outdir = "/hadoop/cms/store/user/jguiang/rare-higgs/{0}_{1}".format(tag, stlTag)
if os.path.isdir("/nfs-7/userdata/jguiang/rare-higgs/stl/{}".format(stlTag)):
    print("Found STL files:")
    print(os.listdir("/nfs-7/userdata/jguiang/rare-higgs/stl/{}".format(stlTag)))
else:
    print("Directory /nfs-7/userdata/jguiang/rare-higgs/stl/{} does not exist.".format(stlTag))

if os.path.isdir("/hadoop/cms/store/user/bemarsh/LGAD/traj_inputs_for_jonathan/{0}/".format(tag)):
    nFiles = len(os.listdir("/hadoop/cms/store/user/bemarsh/LGAD/traj_inputs_for_jonathan/{0}/".format(tag)))
    print("Found {} trajectory files".format(nFiles))
else:
    print("/hadoop/cms/store/user/bemarsh/LGAD/traj_inputs_for_jonathan/{0}/".format(tag))

fout = open("configs/config_{0}_{1}.cmd".format(tag, stlTag), 'w')

username = os.environ["USER"]
x509file = subprocess.check_output(["find","/tmp/", "-maxdepth", "1", "-type", "f", "-user", username, "-regex", "^.*x509.*$"])
os.system("mkdir -p /data/tmp/{0}/condor_job_logs/rare-higgs/".format(username))
os.system("mkdir -p /data/tmp/{0}/condor_submit_logs/rare-higgs/".format(username))

fout.write("""
universe=Vanilla
when_to_transfer_output = ON_EXIT
transfer_input_files=wrapper.sh, package.tar.xz
+DESIRED_Sites="T2_US_UCSD"
+remote_DESIRED_Sites="T2_US_UCSD"
+Owner = undefined
log=/data/tmp/{0}/condor_submit_logs/rare-higgs/condor_12_01_16.log
output=/data/tmp/{0}/condor_job_logs/rare-higgs/1e.$(Cluster).$(Process).out
error =/data/tmp/{0}/condor_job_logs/rare-higgs/1e.$(Cluster).$(Process).err
notification=Never
x509userproxy={1}

""".format(username, x509file))

files = glob.glob("/hadoop/cms/store/user/bemarsh/LGAD/traj_inputs_for_jonathan/{0}/*.txt".format(tag))
jobid = 0
for f in files:
    jobid = (f.split("_")[-1]).split(".")[0]
    fout.write("executable=wrapper.sh\n")
    fout.write("transfer_executable=True\n")
    fout.write("arguments={0} {1} {2} {3}\n".format(jobid, tag, stlTag, outdir))
    fout.write("queue\n\n")


fout.close()
