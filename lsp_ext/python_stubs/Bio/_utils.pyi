from typing import Any, Callable, Protocol, TypeVar

F = TypeVar('F', bound=Callable[..., object])

class _FunctionWithPrevious(Protocol[F]):
    previous: int | None
    __call__: F

def function_with_previous(func: F) -> _FunctionWithPrevious[F]: ...
def find_test_dir(start_dir: str | None = None) -> str: ...
def run_doctest(target_dir: str | None = None, *args: Any, **kwargs: Any) -> None: ...
