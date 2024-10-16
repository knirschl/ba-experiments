import os
import glob
import shutil

def get_metrics(metrics_dir, metric_fname):
  dico = {}
  try:
    with open(os.path.join(metrics_dir, metric_fname)) as writer:
      for line in writer.readlines():
        split = line.split(" : ")
        dico[split[0].lower()] = split[1].replace("\n", "")
  except:
    return None
  return dico

def save_metrics(metrics_dir, dico, metric_fname):
  try:
    os.makedirs(metrics_dir)
  except:
    pass
  with open(os.path.join(metrics_dir, metric_fname), "w") as writer:
    for key, value in sorted(dico.items(), key=lambda x: float(x[1])):
      writer.write(key.lower() + " : " + str(value) + "\n")


def merge():
    metrics_path = "/home/fili/Desktop/2023/BA/code/output/benchmark_results/metrics"
    merge_dir = "-merged"
    metrics = os.listdir(metrics_path)
    metrics = set(z for i in range(len(metrics)) for m in metrics[i+1:] if (z := metrics[i].split('-')[0]) in m) # if starts with same
    for m in metrics:
        try:
            os.mkdir(os.path.join(metrics_path, m + merge_dir))
        except Exception as exc:
            #print(exc)
            pass
        versions = [v for v in glob.glob(os.path.join(metrics_path, m + "*")) if not "-merged" in v]
        merge = []
        for i in range(len(versions)):
            for f1 in os.listdir(versions[i]):
                for v in versions[i+1:]:
                    if f1 not in os.listdir(v):
                        try:
                            shutil.copy(os.path.join(versions[i], f1), os.path.join(metrics_path, m + merge_dir, f1))
                        except Exception as exc:
                            #print(exc)
                            pass
                    else:
                        merge.append(f1)
        dicos = {}
        for f in merge:
            if not f in dicos:
                dicos[f] = []
            for v in versions:
                if ((d := get_metrics(v, f)) != None):
                    dicos[f].append(d)
        #print(dicos)
        merged_dicos = {}
        for f in dicos:
            if f not in merged_dicos:
                merged_dicos[f] = {}
            for dico in dicos[f]:
                for k in dico:
                    if k not in merged_dicos[f]:
                        merged_dicos[f][k] = []
                    merged_dicos[f][k].append(dico[k])
        for f in merged_dicos:
            pops = []
            for k in merged_dicos[f]:
                if k in ["spearfish+fm_full"]: # skip list
                    pops.append(k)
                    continue
                merged_dicos[f][k] = sum([float(e) for e in merged_dicos[f][k]]) / len(merged_dicos[f][k])
            for p in pops:
                merged_dicos[f].pop(p)
            save_metrics(os.path.join(metrics_path, m + merge_dir), merged_dicos[f], f)


def write_over():
    metrics_path = "/home/fili/Desktop/2023/BA/code/output/benchmark_results/metrics"
    merge_dir = "-merged"
    src = os.path.join(metrics_path, "3x80_p")
    src_dist = os.path.join(metrics_path, "1x10_p")
    dest = src + merge_dir
    try:
        os.mkdir(dest)
    except Exception as exc:
        print(exc)
        pass
    for f in os.listdir(src):
        if not "rf_dist" in f:
            try:
                shutil.copy(os.path.join(src, f), os.path.join(dest, f))
                continue
            except Exception as exc:
                print(exc)
                pass
        dico_src = get_metrics(src, f)
        dico_src_dist = get_metrics(src_dist, f)
        dico = {}
        for k in dico_src:
            if "spearfish.p_" in k:
                dico[k] = dico_src[k]
            else:
                dico[k] = dico_src_dist[k]
        save_metrics(dest, dico, f)

#write_over()