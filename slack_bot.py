import requests
from pathlib import Path

def send_message(message: str, channel: str="C0B0TD3LD3N"):
    DATA_DIR = Path(__file__).parent 
    BOT_TOKEN: str = (DATA_DIR / "slack_api.key").read_text().strip()
    headers: dict = {
    "Authorization": "Bearer " + BOT_TOKEN,
    "Content-Type": "application/json",
    }
    
    json_data: dict = { 
    "channel": "",
    "text": "",
    }
    json_data["text"] = message  # set the text field to the given message
    json_data["channel"] = (
        channel  # set the channel field to the given channel
    )

    requests.post(
        "https://slack.com/api/chat.postMessage", headers=headers, json=json_data
    ) # Sends the message to slack servers

# if __name__ == "__main__":
        # PAPER = "Sina"
        # dataset_name = "PoreMaker"
        # targ = "permeability (GPU)"

        # send_message(f"""
        # Finished training
        # PAPER: {PAPER}
        # Dataset: {dataset_name}
        # Target: {targ}
        # Regressor: GPytorchMAP
        # """)