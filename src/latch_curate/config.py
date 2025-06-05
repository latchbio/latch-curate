from pathlib import Path

from latch_curate.constants import latch_curate_constants

class UserConfig:

    def __init__(self):
        self._root = Path.home().resolve() / ".latch"
        self._token_path = None

    @property
    def root(self) -> Path:
        if not self._root.exists():
            self._root.mkdir(parents=True)
        return self._root

    @property
    def package_version_cache_location(self) -> Path:
        return self.root / latch_curate_constants.pkg_version_cache_path

    @property
    def openai_api_key(self) -> Path:
        key_file = self.root / latch_curate_constants.openai_api_key_path
        assert key_file.exists()
        return key_file.read_text().strip()

user_config = UserConfig()
