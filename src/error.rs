//! Error types for BWA-MEM operations.

use thiserror::Error;

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