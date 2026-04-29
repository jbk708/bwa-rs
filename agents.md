# BWA-MEM Coding Standards

## Principles
1. Zero-cost abstractions, no unnecessary heap allocation
2. Memory-efficient: 2-bit encoding, memory mapping for 3GB genomes
3. Minimal dependencies - implement from scratch
4. No panics in library code - use `Result<T, BwaError>`
5. Code speaks for itself - no inline comments explaining *what*

## Workflow (One Ticket Per PR)

1. **Branch** - Create feature branch: `git checkout -b t{N}-short-description`
2. **Plan** - Read ticket, outline types/functions/tests, get approval
3. **Tests** - Write failing tests first, define expected behavior
4. **Implement** - Minimum code to pass tests
5. **Simplify** - `cargo clippy --fix`, no warnings, `cargo test` passes
6. **PR** - Create PR: `gh pr create --title "T{N}: {Title}" --body "..."`

## Formats

### Ticket
```markdown
### T{N}: {Title}
**Description:** One sentence

**Deliverables:**
- [ ] Item 1

**Dependencies:** T{n}
```

### Commit
```
T{N}: {Verb} {what changed}
```

### PR Title
```
T{N}: {Ticket title}
```

### PR Body
Use single quotes around 'gh pr create --body' content to avoid shell interpolation issues with double quotes.

## Naming
- Types: `CamelCase` - `Sequence`, `FMIndex`
- Functions/variables: `snake_case` - `find_mems`, `insert_size`
- Constants: `SCREAMING_SNAKE_CASE`

## Conventions
- Use slices (`&[u8]`) over `Vec` for borrowing
- Prefer `Result` for fallible operations
- Return iterators for streaming large data
- 2-bit encoding: A=0, C=1, G=2, T=3, N=4
- Positions: 0-indexed in code, 1-indexed in SAM

## Error Type
```rust
#[derive(Error, Debug)]
pub enum BwaError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    Parse(String),
    #[error("Index error: {0}")]
    Index(String),
    #[error("Alignment error: {0}")]
    Alignment(String),
}
```

## Review Checklist
- [ ] `cargo test` passes
- [ ] `cargo clippy` clean
- [ ] No panics, all fallibles return `Result`
- [ ] No unnecessary clones or comments
- [ ] Tests cover edge cases