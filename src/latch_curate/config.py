import json
from pathlib import Path

from latch_curate.constants import latch_curate_constants

class EmailConfig():
    def __init__(self, config_path: Path | str | None = None):
        default_path = Path.home() / ".latch" / "latch-curate" / "email-info.json"
        self._config_file = Path(config_path) if config_path else default_path

        config_data: dict = {}
        if self._config_file.is_file():
            try:
                with self._config_file.open("r") as f:
                    config_data = json.load(f)
            except json.JSONDecodeError:
                raise RuntimeError(f"Invalid JSON in {self._config_file}")

        try:
            self.smtp_host     = config_data["smtp_host"]
            self.smtp_port     = config_data["smtp_port"]
            self.smtp_user     = config_data["smtp_user"]
            self.smtp_password = config_data["smtp_password"]
            self.sender_addr   = config_data["sender_addr"]
            self.starttls      = config_data["starttls"]
            self.timeout       = config_data["timeout"]
        except KeyError:
            # todo(kenny)
            print("No email config detected: limits publish functionality")
            # raise RuntimeError(f"Malformed {self._config_file}")

email_config = EmailConfig()

class UserConfig:

    def __init__(self):
        self._root = Path.home().resolve() / ".latch"

    @property
    def root(self) -> Path:
        if not self._root.exists():
            self._root.mkdir(parents=True)
        return self._root

    @property
    def token(self) -> Path:
        token_path = self.root / "token"
        if not token_path.exists():
            raise ValueError("SDK token does not exist. Ensure you are in a properly configured latch pod.")
        return token_path.read_text()

    @property
    def package_version_cache_location(self) -> Path:
        return self.root / latch_curate_constants.pkg_version_cache_path

    @property
    def openai_api_key(self) -> Path:
        key_file = self.root / latch_curate_constants.openai_api_key_path
        assert key_file.exists()
        return key_file.read_text().strip()

user_config = UserConfig()
