try:
    import GPUtil
except ImportError(GPUtil):
    raise ModuleNotFoundError('please install GPUtil to use this SOP with "pip install GPUtil -y"')

import GPUtil
import os
import argparse


def getpid():
    mypid = os.getpid()
    pidFile = open(".mypid", "w")
    pidFile.write(str(mypid))
    pidFile.close()


def getGPUids(_excluded):
    GPUs = GPUtil.getGPUs()
    gpu_ids = []
    for availableGPU in GPUs:
        gpu_ids.append(availableGPU.id)
    if _excluded:
        for excluded_GPU in _excluded:
            gpu_ids.remove(int(excluded_GPU))
    if len(gpu_ids) != 0:
        print("Available GPUS: ", gpu_ids)
        return gpu_ids
    else:
        print("Please leave at least one GPU to run mwSuMD and run again.")
        exit()


def createBatches(b_replicas, b_total_gpu_ids):
    quotient, rest = divmod(b_replicas, len(b_total_gpu_ids))
    b_result = quotient * b_total_gpu_ids + b_total_gpu_ids[:rest]
    b_batches = [b_result[b_i:b_i + len(b_total_gpu_ids)] for b_i in range(0, len(b_result), len(b_total_gpu_ids))]
    return b_batches, b_total_gpu_ids


def DescribeSettings(arguments):
    print("Settings:")
    for arg, value in arguments.items():
        print(arg, " ==> ", value)
    print("\n")


class ParseMembraneRestraints(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        membrane_list = []
        if values:
            for pair in values:
                parts = pair.split(",")
                if len(parts) == 2:
                    membrane_list.append(parts)
                else:
                    raise argparse.ArgumentTypeError(f"{pair} is not in the correct RESNAME,ATOM format.")
        setattr(namespace, self.dest, membrane_list)
