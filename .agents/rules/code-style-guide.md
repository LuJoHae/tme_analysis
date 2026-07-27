---
trigger: always_on
---

# Comprehensive Functional Python Coding Style Rules

This document outlines the strict coding standards for functional programming in Python. We enforce a declarative, pure, and immutable approach to software design to maximize predictability, testability, and composability.

## 1. Core Functional Principles
- **Pure Functions**: Functions must map inputs to outputs without any side effects. They should not mutate global state, write to files (unless explicitly handled via IO boundaries), or modify incoming arguments.
- **Declarative & Expression-Oriented Code**: Describe *what* to do rather than *how* to do it. Treat the codebase as expression-oriented: prefer ternary operators (`x if cond else y`) or `match` over `if/else` statements. The `global`, `nonlocal`, and `del` keywords are strictly forbidden.
- **First-Class Functions**: Treat functions as data. Use higher-order functions to compose logic.
- **Functional Core, Imperative Shell**: Architect the application so the "Core" consists solely of pure functions. Push all side effects (database I/O, network requests, logging) to the "Shell" (the outer boundaries of the application).

## 2. Architecture & Design Patterns
- **Typeclasses over Inheritance**: Do not use Object-Oriented inheritance (`class A(B):`). Strictly separate data (frozen models) from behavior (functions). Implement polymorphic behavior using **Typeclasses** via `typing.Protocol` (duck typing).
- **Dependency Injection**: Functions should never instantiate their own side-effecting dependencies (like database clients). Pass dependencies as arguments to maintain function purity and ease of testing.
- **Structural Pattern Matching**: Use Python 3.10+ `match...case` statements to declaratively unpack data structures and route logic.
- **Exhaustive Pattern Matching**: All `match...case` blocks must be exhaustive. Use `typing.assert_never(value)` in the default `case _:` block so static analysis (`mypy`) fails if a state or enum case is unhandled.
- **Lazy Evaluation**: Prefer generator expressions `(...)` and the `yield` keyword over list comprehensions `[...]` when processing large datasets.

## 3. Error Handling & Railway Oriented Programming
- **No Exceptions for Control Flow**: Functions should **not** throw exceptions (`raise`) for expected business logic failures. 
- **Return Errors (Result Monad)**: Return errors as values using a `Result` pattern. 
- **Monadic Do-Notation & Pipelining**: Avoid manually unwrapping `Result` types. Chain operations using `.map()` and `.bind()`. Alternatively, use **Do-Notation** (via the `@do` decorator from `returns`) to write sequential monadic code that flattens callback hell while remaining mathematically pure.

## 4. Strict Immutability & Data Modeling
- **Pydantic for Validation**: Use `pydantic` models configured with `frozen=True` (or `@dataclass(frozen=True)`) for data modeling. 
- **No In-Place Mutation**: Never modify data in place.
- **Lenses for Deep Updates**: When updating deeply nested immutable data structures, use a **Lenses** library (e.g., `lenses` or `optics`) to perform the update cleanly instead of writing verbose `.model_copy(update=...)` calls.
- **Ban `None`**: Avoid the use of `None` and `Optional[T]` entirely. Instead, use a `Maybe[T]` or `Option[T]` monad to force explicit handling of missing values.

## 5. Declarative Data Processing & Visualization
- **Dataframes**: Use **Polars** for all dataframe operations instead of `pandas`. Polars provides a declarative, lazy-evaluation API that prevents inplace mutations.
- **Plotting Library**: Use **Altair** (based on Vega/Vega-Lite) for all plotting. Do not use Matplotlib or Seaborn.
- **Plot Export**: Plots should **not** be saved as HTML. Instead, export them as **SVG** (Scalable Vector Graphics). SVG is the recommended format because it preserves high-quality, resolution-independent vector graphics that are ideal for publications and detailed analysis. (Note: exporting static images from Altair requires `vl-convert-python`).

## 6. Testing & Quality Assurance
- **Property-Based Testing**: Use `hypothesis` to write property-based tests. Define invariants and let the framework generate randomized inputs to aggressively stress-test your pure functions.
- **Strict Static Typing**: The codebase must pass static analysis using **`mypy --strict`**.
- **Specialized FP Linters**: The CI/CD pipeline must enforce functional purity using specialized linters. We strongly recommend configuring **`wemake-python-styleguide`**, which strictly bans variable shadowing, side-effects, and mutable state by default.

## Example
```python
from typing import Generic, TypeVar, Protocol, assert_never
from pydantic import BaseModel, ConfigDict
from returns.result import Result, Success, Failure
from returns.maybe import Maybe, Some, Nothing
from returns.pipeline import is_successful
from returns.methods import bind
from returns.pointfree import map_
from returns.decorators import do
import altair as alt
import polars as pl

# 1. Immutable Data Model (No None allowed)
class UserData(BaseModel):
    model_config = ConfigDict(frozen=True)
    id: int
    score: float
    category: Maybe[str]

# 2. Typeclass (Protocol) instead of Inheritance
class DatabaseClient(Protocol):
    def fetch_users(self) -> Result[list[dict], str]: ...

# 3. Pure Function & Expression-Oriented Pattern Matching
def categorize_score(score: float) -> str:
    match score:
        case s if s > 90:
            return "A"
        case s if s > 80:
            return "B"
        case s if s <= 80:
            return "C"
        case _ as unreachable:
            assert_never(unreachable)

# 4. Pure Data Processing (Polars)
def process_user_data(data: tuple[UserData, ...]) -> Result[pl.DataFrame, str]:
    return Failure("No data.") if not data else Success(
        pl.DataFrame([
            {"id": u.id, "score": u.score, "category": u.category.value_or("Unknown")} 
            for u in data
        ]).lazy()
        .with_columns(
            pl.col("score").map_elements(categorize_score, return_dtype=pl.String).alias("grade")
        ).collect()
    )

# 5. Declarative Plotting (Altair)
def create_score_chart(df: pl.DataFrame) -> Result[alt.Chart, str]:
    return Failure("Empty.") if df.is_empty() else Success(
        alt.Chart(df).mark_bar().encode(
            x='grade:N',
            y='count():Q',
            color='grade:N'
        ).properties(title="Grades")
    )

# 6. Monadic Do-Notation Pipeline (Imperative Shell)
@do(Result[alt.Chart, str])
def run_pipeline(db_client: DatabaseClient) -> alt.Chart:
    """Uses Do-Notation to unwrap Success values implicitly; short-circuits on Failure."""
    # Yield unwraps the Result if Success, or returns early if Failure
    raw_data = yield db_client.fetch_users()
    
    # Parse data
    users = tuple(
        UserData(id=r["id"], score=r["score"], category=Maybe.from_optional(r.get("category"))) 
        for r in raw_data
    )
    
    # Process and Plot
    df = yield process_user_data(users)
    chart = yield create_score_chart(df)
    
    return chart

def main(db_client: DatabaseClient) -> None:
    match run_pipeline(db_client):
        case Failure(err):
            print(f"Pipeline failed: {err}")
        case Success(chart):
            print("Chart generated successfully.")
```
