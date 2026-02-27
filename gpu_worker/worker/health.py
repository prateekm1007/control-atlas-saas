"""
Toscanini Phase B2 — Worker Health Monitor
"""
import os, time

REDIS_URL = os.environ.get("REDIS_URL", "redis://redis:6379/0")


def check_redis() -> dict:
    try:
        import redis
        r = redis.from_url(REDIS_URL, socket_connect_timeout=3)
        t = time.time()
        r.ping()
        info = r.info("server")
        return {
            "status":          "ok",
            "latency_ms":      round((time.time() - t) * 1000, 2),
            "redis_version":   info.get("redis_version", "unknown"),
            "used_memory":     info.get("used_memory_human", "unknown")
        }
    except Exception as e:
        return {"status": "error", "error": str(e)}


def check_gpu() -> dict:
    try:
        import subprocess
        r = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total,memory.free",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=5
        )
        if r.returncode == 0:
            gpus = []
            for line in r.stdout.strip().split("\n"):
                p = [x.strip() for x in line.split(",")]
                if len(p) == 3:
                    gpus.append({"name": p[0],
                                 "memory_total_mb": int(p[1]),
                                 "memory_free_mb":  int(p[2])})
            return {"status": "ok", "gpu_count": len(gpus), "gpus": gpus}
        return {"status": "no_gpu"}
    except Exception as e:
        return {"status": "no_gpu", "error": str(e)}


def check_openmm() -> dict:
    try:
        import openmm as mm
        platforms = [mm.Platform.getPlatform(i).getName()
                     for i in range(mm.Platform.getNumPlatforms())]
        return {
            "status":         "ok",
            "version":        mm.__version__,
            "platforms":      platforms,
            "cuda_available": "CUDA" in platforms
        }
    except Exception as e:
        return {"status": "error", "error": str(e)}


def full_health_check() -> dict:
    return {
        "redis":        check_redis(),
        "gpu":          check_gpu(),
        "openmm":       check_openmm(),
        "worker_ready": True
    }


if __name__ == "__main__":
    import json
    print(json.dumps(full_health_check(), indent=2))
