use super::common::Status;


pub type AlignmentResult = Result<(i32, Status), AlignmentError>;

#[derive(thiserror::Error, Debug)] 
pub enum AlignmentError {
    #[error("{0}")]
    InvalidRangeError(String),
    #[error("{0}")]
    QueryRangeError(String),
    #[error("{0}")]
    ReferenceRangeError(String),
    #[error("{0}")]
    InvalidAlignmentError(String),
}