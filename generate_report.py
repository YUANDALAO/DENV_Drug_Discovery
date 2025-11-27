import pandas as pd
from tensorboard.backend.event_processing import event_accumulator
import glob
import os
import matplotlib.pyplot as plt
from datetime import datetime

# 读取TensorBoard数据
log_dir = "experiments/runs/run15/tensorboard_0"
event_files = glob.glob(os.path.join(log_dir, "events.out.*"))

if event_files:
    latest = max(event_files, key=os.path.getmtime)
    ea = event_accumulator.EventAccumulator(latest)
    ea.Reload()
    
    # 生成HTML报告
    html = f"""
    <html>
    <head><title>DENV Training Report - {datetime.now().strftime('%Y-%m-%d')}</title></head>
    <body>
    <h1>REINVENT4 Training Report</h1>
    <h2>Run: run15</h2>
    <h3>Generated: {datetime.now()}</h3>
    """
    
    # 添加所有指标
    for tag in ea.Tags()['scalars']:
        events = ea.Scalars(tag)
        if events:
            latest_event = events[-1]
            html += f"<p><b>{tag}</b>: {latest_event.value:.4f} (Step {latest_event.step})</p>"
    
    # 添加候选分子
    if os.path.exists("experiments/runs/run15/candidates_gold.csv"):
        df = pd.read_csv("experiments/runs/run15/candidates_gold.csv")
        html += "<h2>Top Candidates</h2>"
        html += df.head(10).to_html()
    
    html += "</body></html>"
    
    with open("training_report.html", "w") as f:
        f.write(html)
    
    print("✓ 报告已生成: training_report.html")
