import os
import time
import threading
import arcpy
from arcpy.sa import Aggregate, Raster

arcpy.CheckOutExtension("Spatial")
arcpy.env.overwriteOutput = True

DEM_30M = r"Z:\download\Mamad\Europe\SubCatch_from_DEM30m\gedtm_rf_m_30m_s_2000_22_3035_v2025.tif"
OUT_ROOT = r"E:\pythonProject\Europe"
OUT_DIR = os.path.join(OUT_ROOT, "PREP_90M")
os.makedirs(OUT_DIR, exist_ok=True)

DEM_90M = os.path.join(OUT_DIR, "dem_90m_3035.tif")

AGG_FACTOR = 3
AGG_TYPE = "MINIMUM"

SCRATCH = os.path.join(OUT_DIR, "_scratch")
os.makedirs(SCRATCH, exist_ok=True)
arcpy.env.scratchWorkspace = SCRATCH


class StepTimer:
    def __init__(self):
        self._t0 = None

    def start(self, label: str):
        self._t0 = time.perf_counter()
        print(f"[START] {label}")

    def stop(self, label: str):
        dt = time.perf_counter() - self._t0
        print(f"[DONE ] {label} in {dt/60:.2f} min ({dt:.1f} s)")
        return dt


def _file_growth_monitor(path: str, stop_event: threading.Event, interval=10):
    """
    Prints output file size while a tool is running.
    Useful as a rough 'is it doing something?' indicator and to infer ETA after it finishes.
    """
    last_size = -1
    t0 = time.perf_counter()
    while not stop_event.is_set():
        try:
            size = os.path.getsize(path) if os.path.exists(path) else 0
        except OSError:
            size = 0
        if size != last_size:
            elapsed = time.perf_counter() - t0
            print(f"  [WRITE] {os.path.basename(path)} size={size/1e9:.2f} GB elapsed={elapsed/60:.1f} min")
            last_size = size
        stop_event.wait(interval)


def main():
    if not os.path.exists(DEM_30M):
        raise FileNotFoundError(f"DEM_30M not found: {DEM_30M}")

    timer = StepTimer()

    # Optional: if you have a typical machine and only run once, an "estimate" can be based on input size.
    in_gb = os.path.getsize(DEM_30M) / 1e9
    print(f"Input size: {in_gb:.2f} GB (30m). Output will be much smaller (~1/9 pixels).")

    # Step: Aggregate -> 90m
    label = f"Aggregate {AGG_FACTOR}x ({AGG_TYPE}) to 90m"
    timer.start(label)

    # Monitor output growth during write (approximate progress signal)
    stop_evt = threading.Event()
    mon = threading.Thread(target=_file_growth_monitor, args=(DEM_90M, stop_evt), daemon=True)
    mon.start()

    try:
        dem90 = Aggregate(Raster(DEM_30M), AGG_FACTOR, AGG_TYPE, "EXPAND", "DATA")
        dem90.save(DEM_90M)
    finally:
        stop_evt.set()
        mon.join(timeout=2)

    dt_agg = timer.stop(label)

    # Step: CalculateStatistics
    label2 = "CalculateStatistics"
    timer.start(label2)
    arcpy.management.CalculateStatistics(DEM_90M)
    dt_stats = timer.stop(label2)

    # Report output raster properties
    r = Raster(DEM_90M)
    print("Output:", DEM_90M)
    print("Cellsize:", r.meanCellWidth, r.meanCellHeight)
    print("CRS:", r.spatialReference.factoryCode, r.spatialReference.name)
    print("Rows/Cols:", r.height, r.width)

    # Simple “next run” estimate based on this run’s timings:
    # (This is the only reliable ETA you can get: learn from your machine.)
    print("\n--- Timing summary ---")
    print(f"Aggregate: {dt_agg/60:.2f} min")
    print(f"Stats:     {dt_stats/60:.2f} min")
    print(f"Total:     {(dt_agg+dt_stats)/60:.2f} min")

if __name__ == "__main__":
    main()
