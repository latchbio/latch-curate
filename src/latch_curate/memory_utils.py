import os
import resource
from pathlib import Path
from typing import Optional


def get_memory_limit() -> int:
    cgroup_limit = _get_cgroup_memory_limit()
    if cgroup_limit is not None:
        return cgroup_limit

    root_limit = _read_memory_max("/sys/fs/cgroup/memory.max")
    if root_limit is not None:
        return root_limit

    return _get_system_memory()


def _get_cgroup_memory_limit() -> Optional[int]:
    try:
        cgroup_path = Path("/proc/self/cgroup")
        if not cgroup_path.exists():
            return None

        cgroup_content = cgroup_path.read_text().strip()

        for line in cgroup_content.splitlines():
            if line.startswith("0::"):
                relative_path = line.split("::", 1)[1]
                if relative_path == "/":
                    relative_path = ""
                break
        else:
            return None

        base_path = Path("/sys/fs/cgroup")
        min_limit = None

        if relative_path:
            path_parts = relative_path.strip("/").split("/")

            for i in range(len(path_parts) + 1):
                if i == 0:
                    check_path = base_path
                else:
                    check_path = base_path / "/".join(path_parts[:i])

                limit = _read_memory_max(check_path / "memory.max")
                if limit is not None:
                    if min_limit is None or limit < min_limit:
                        min_limit = limit
        else:
            min_limit = _read_memory_max(base_path / "memory.max")

        return min_limit

    except Exception:
        return None


def _read_memory_max(path: Path) -> Optional[int]:
    try:
        if isinstance(path, str):
            path = Path(path)

        if not path.exists():
            return None

        content = path.read_text().strip()

        if content == "max":
            return None

        return int(content)

    except (IOError, ValueError):
        return None


def _get_system_memory() -> int:
    try:
        pages = os.sysconf("SC_PHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        return pages * page_size
    except (ValueError, AttributeError):
        pass

    try:
        meminfo = Path("/proc/meminfo")
        if meminfo.exists():
            for line in meminfo.read_text().splitlines():
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    return kb * 1024
    except Exception:
        pass

    try:
        soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_AS)
        if hard_limit != resource.RLIM_INFINITY:
            return hard_limit
    except Exception:
        pass

    return 8 * 1024 * 1024 * 1024


def format_memory(bytes_value: int) -> str:
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if bytes_value < 1024:
            return f"{bytes_value:.2f} {unit}"
        bytes_value /= 1024
    return f"{bytes_value:.2f} PB"


if __name__ == "__main__":
    limit = get_memory_limit()
    print(f"Detected memory limit: {format_memory(limit)}")

    cgroup_limit = _get_cgroup_memory_limit()
    if cgroup_limit:
        print(f"  - Cgroup limit: {format_memory(cgroup_limit)}")

    root_limit = _read_memory_max("/sys/fs/cgroup/memory.max")
    if root_limit:
        print(f"  - Root cgroup limit: {format_memory(root_limit)}")

    system_mem = _get_system_memory()
    print(f"  - System memory: {format_memory(system_mem)}")
