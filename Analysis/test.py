from pathlib import Path
import pandas as pd

# Define the weekly strength plan
data = {
    "Day": ["Monday", "Tuesday", "Thursday", "Saturday", "Sunday"],
    "Focus": [
        "Upper Body #1 (Push-Up Progression – Home)",
        "Leg Day (Post-Run)",
        "Upper Body #2 (Push-Up Progression – Home)",
        "Gym Strength (Full Upper Body)",
        "Recovery / Optional Mobility"
    ],
    "Duration": ["20 min", "20–30 min", "20 min", "45 min", "Optional"],
    "Details": [
        "Incline + Eccentric + Kneeling Push-Ups, Dumbbell Rows/Press, Plank",
        "Glute Bridges, Split Squats, Wall Sits, Band Walks, Calf Raises",
        "Repeat Mon with added reps or lower incline",
        "Incline/Full Push-Ups, Bench/Chest Press, Rows, Overhead Press, Core",
        "Mobility, light yoga, or complete rest"
    ]
}

# Create DataFrame
df = pd.DataFrame(data)

# Save as Excel and CSV for Notion compatibility
excel_path = "/mnt/data/Weekly_Strength_Training_Plan_April2025.xlsx"
csv_path = "/mnt/data/Weekly_Strength_Training_Plan_April2025.csv"
df.to_excel(excel_path, index=False)
df.to_csv(csv_path, index=False)

import ace_tools as tools; tools.display_dataframe_to_user(name="Weekly Strength Training Plan", dataframe=df)

excel_path, csv_path
