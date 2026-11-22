// The world data file is built by .github/workflows/data.yml and published as
// a GitHub release; it is not kept in git. Download it to data/world.bin
// before building (CI and release workflows do this automatically).
fn main() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/world.bin");
    println!("cargo:rerun-if-changed={}", path.display());
    if !path.exists() {
        panic!(
            "\n\nMissing {}.\nDownload it with (N = the number in tools/data/VERSION):\n  curl -L -o data/world.bin \
             https://github.com/djeyprott/sim3/releases/download/world-data-vN/world.bin\n\n",
            path.display()
        );
    }
}
