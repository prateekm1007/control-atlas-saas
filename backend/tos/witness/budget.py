def enforce_narrative_budget(text, limit=1000):
    words = text.split()
    return " ".join(words[:limit]) + (" [TRUNCATED]" if len(words) > limit else "")
