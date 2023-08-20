import sys
import os
import shutil
import subprocess
import paths

def str_2(ll):
  return "{0:.2f}".format(ll)

def str_4(ll):
  return "{0:.4f}".format(ll)

def printFlush(msg):
  print(msg)
  sys.stdout.flush()

def relative_symlink(src, dest):
  relative_path = os.path.relpath(src, os.path.dirname(dest))
  tmp = dest + ".sym"
  os.symlink(relative_path,  tmp)
  shutil.move(tmp, dest)

def reset_dir(dir_name):
  shutil.rmtree(dir_name, True)
  os.makedirs(dir_name)

def create_result_dir(suffix, additional_args = []):
  base = os.path.join(results_root, suffix) 
  for arg in additional_args:
    base += "_" + arg
  base += "_"
  result_dir = ""
  for i in range(0, 10000):
    result_dir = base + str(i)
    if (not os.path.isdir(result_dir)):
      os.makedirs(result_dir)
      #open(historic, "a+").write("Results directory: " + result_dir + "\n")
      #print("Results directory: " + result_dir)
      return os.path.abspath(result_dir)

def checkAndDelete(arg, arguments):
  if (arg in arguments):
    arguments.remove(arg)
    return True
  return False

def getAndDelete(arg, arguments, default_value):
  print ("looking for " + arg + " in " + str(arguments))
  if (arg in arguments):
    index = arguments.index(arg)
    res = arguments[index + 1]
    del arguments[index + 1]
    del arguments[index]
    return res
  else:
    return default_value

def getArg(arg, arguments, default_value):
  print ("looking for " + arg + " in " + str(arguments))
  if (arg in arguments):
    index = arguments.index(arg)
    res = arguments[index + 1]
    return res
  else:
    return default_value

def submit_normal(submit_file_path, submit_id, command, log_cout):
    commands_list = command.split("\n")
    logfile = os.path.join(os.path.dirname(submit_file_path), submit_id + "_logs.out")
    for subcommand in commands_list:
      if (log_cout):
        subprocess.check_call(subcommand, shell=True)
      else:
        subprocess.check_call(subcommand + " &>> " + logfile , shell=True)

def submit_haswell(submit_file_path, submit_id, command, threads, debug):
  threads = int(threads)
  nodes = str((int(threads) - 1) // 16 + 1)
  logfile = os.path.join(os.path.dirname(submit_file_path), submit_id + "_logs.out")
  with open(submit_file_path, "w") as f:
    f.write("#!/bin/bash\n")
    f.write("#SBATCH -o " + logfile + "\n")
    f.write("#SBATCH -B 2:8:1\n")
    f.write("#SBATCH -N " + str(nodes) + "\n")
    f.write("#SBATCH -n " + str(threads) + "\n")
    f.write("#SBATCH --threads-per-core=1\n")
    f.write("#SBATCH --cpus-per-task=1\n")
    f.write("#SBATCH --hint=compute_bound\n")
    if (debug):
      f.write("#SBATCH -t 2:00:00\n")
    else:
      f.write("#SBATCH -t 24:00:00\n")

    f.write("\n")
    f.write(command)
  command = []
  command.append("sbatch")
  if (debug):
    command.append("--qos=debug")
  command.append("-s")
  command.append(submit_file_path)
  out = open(paths.historic, "a+")
  subprocess.check_call(command, stdout = out)
  out.write("Output in " + logfile + "\n")
  print(open(paths.historic).readlines()[-1][:-1])
  out.write("\n")

def submit_cascade(submit_file_path, submit_id, command, threads, debug):
  threads = int(threads)
  nodes = str((int(threads) - 1) // 20 + 1)
  logfile = os.path.join(os.path.dirname(submit_file_path), submit_id + "_logs.out")
  with open(submit_file_path, "w") as f:
    f.write("#!/bin/bash\n")
    f.write("#SBATCH -o " + logfile + "\n")
    #f.write("#SBATCH -B 2:8:1\n")
    f.write("#SBATCH -N " + str(nodes) + "\n")
    f.write("#SBATCH -n " + str(threads) + "\n")
    f.write("#SBATCH --threads-per-core=1\n")
    f.write("#SBATCH --cpus-per-task=20\n")
    f.write("#SBATCH --hint=compute_bound\n")
    if (debug):
      f.write("#SBATCH -t 2:00:00\n")
    else:
      f.write("#SBATCH -t 24:00:00\n")

    f.write("\n")
    f.write(command)
  command = []
  command.append("sbatch")
  if (debug):
    command.append("--qos=debug")
  command.append("-s")
  command.append(submit_file_path)
  out = open(paths.historic, "a+")
  subprocess.check_call(command, stdout = out)
  out.write("Output in " + logfile + "\n")
  print(open(paths.historic).readlines()[-1][:-1])
  out.write("\n")

def submit(submit_file_path, command, threads, cluster):
  submit_id = os.path.basename(submit_file_path).replace("run_", '').replace(".sh", '')
  if (cluster == "normal"):
    submit_normal(submit_file_path, submit_id, command, False)
  elif (cluster == "normald"):
    submit_normal(submit_file_path, submit_id, command, True)
  elif (cluster == "haswell"):
    submit_haswell(submit_file_path, submit_id, command, threads, False)
  elif (cluster == "haswelld"):
    submit_haswell(submit_file_path, submit_id, command, threads, True)
  elif (cluster == "cascade"):
    submit_cascade(submit_file_path, submit_id, command, threads, False)
  elif (cluster == "cascaded"):
    submit_cascade(submit_file_path, submit_id, command, threads, True)
  elif (cluster == "magny"):
    submit_magny(submit_file_path, submit_id, command, threads)
  else:
    print("unknown cluster " + cluster)
    sys.exit(1)

def run_with_scheduler(executable, command_file, parallelization, cores, scheduler_output_dir, logname = None):
  command = ""
  out = sys.stdout
  if (logname != None):
    out = open(os.path.join(scheduler_output_dir, logname), "w")

  isMPI = (parallelization == "onecore") or (parallelization == "split")
  if (isMPI):
    command += "mpirun -np " + str(cores) + " "
  command += paths.mpischeduler_exec + " "
  command += "--" + parallelization + "-scheduler "
  command += str(cores) + " "
  command += executable + " "
  command += command_file + " "
  command += scheduler_output_dir 
  print("Running " + command)
  subprocess.check_call(command.split(" "), stdout = out, stderr = out)
