import google.generativeai as genai
import os

# Using the Sovereign Key
KEY = "AIzaSyD1d8ydNS345g_qblFNaSAntBwifk_L-Bg"
genai.configure(api_key=KEY)

print("\n--- 🛡️ TOSCANINI APEX CAPABILITY PROBE ---")
try:
    available = [m.name for m in genai.list_models() if 'generateContent' in m.supported_generation_methods]
    for m in available:
        print(f"✅ AUTHORIZED: {m}")
except Exception as e:
    print(f"🚨 PROBE FAILED: {str(e)}")
