import resource, subprocess, sys, os
p = subprocess.run([sys.executable, "score204.py", "--arms", "equal", "a_self",
                    "--folds", "orion_hek293t", "--out", "work/scores204-rssprobe.json"],
                   cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
r = resource.getrusage(resource.RUSAGE_CHILDREN)
print(f"peak child RSS {r.ru_maxrss/(1<<30):.3f} GiB")
