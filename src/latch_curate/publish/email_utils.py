import smtplib
import ssl
import mimetypes
from email.message import EmailMessage
from pathlib import Path
from dataclasses import dataclass

from latch_curate.config import email_config

@dataclass
class EmailRecipient:
    email: str
    name: str

def check_config():
    missing = [k for k in
               ("smtp_host","smtp_port","smtp_user","smtp_password","sender_addr")
               if getattr(email_config, k) in (None, "")]
    if len(missing) > 0:
        raise ValueError(f"Missing email configs: {missing}")

def attach_files(msg: EmailMessage,
                  attachments: list[Path]):
    for p in attachments:
        p = Path(p)
        ctype, encoding = mimetypes.guess_type(p)
        maintype, subtype = (ctype or "application/octet-stream").split("/", 1)
        with p.open("rb") as fh:
            msg.add_attachment(fh.read(),
                               maintype=maintype,
                               subtype=subtype,
                               filename=p.name)

def send_email_to_authors(subject: str, body: str, recipients: list[EmailRecipient], report_path: Path):
    check_config()
    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"] = email_config.sender_addr
    msg["To"] = ", ".join([x.email for x in recipients])

    print(body)

    msg.set_content(body, subtype="html")

    attach_files(msg, [report_path])

    context = ssl.create_default_context()
    with smtplib.SMTP(email_config.smtp_host, email_config.smtp_port, timeout=email_config.timeout) as s:
        if email_config.starttls:
            s.starttls(context=context)
        s.login(email_config.smtp_user, email_config.smtp_password)
        s.send_message(msg)
    print(f"message sent to {recipients}")
