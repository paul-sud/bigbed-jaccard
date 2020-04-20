use std::error::Error;
use std::io::Write;

use reqwest;
use tempfile::NamedTempFile;

/// Asynchronously download to a temporary file. Using tokio::main here so callers can
/// just treat this function as sync
#[tokio::main]
pub async fn download_to_tempfile(url: &str) -> Result<NamedTempFile, Box<dyn Error>> {
    let mut tmpfile = NamedTempFile::new()?;
    let response = reqwest::get(url).await?;
    tmpfile.write_all(&response.bytes().await?)?;
    Ok(tmpfile)
}
