# 程序所需用到的一些工具
import time
from functools import wraps


def print_runtime(func):
    """打印函数运行时间的装饰器（仅支持同步函数）"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)  # 执行原函数
        end_time = time.perf_counter()
        print(f"函数 {func.__name__} 运行耗时: {end_time - start_time:.6f} 秒")
        return result
    return wrapper


