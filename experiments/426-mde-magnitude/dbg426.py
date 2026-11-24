"""426 debug: run the MDE arm target by target and report the first placement failure's class balance."""
import sys
sys.argv = ["run426.py", "base", "9", "--bulk", "ps", "--mde-c", "0.1"]
import runpy
runpy.run_path("run426.py", run_name="__main__")
