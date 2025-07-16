import os

class DataManager:
    def get_instance(self):
        if ("THE_BIOM_SIGNATURES" in os.environ and "THE_BIOM_EXPRESSIONS" in os.environ and "THE_BIOM_PATHWAYS" in os.environ):
            signatures = os.environ["THE_BIOM_SIGNATURES"].strip()
            expressions = os.environ["THE_BIOM_EXPRESSIONS"].strip()
            pathways = os.environ["THE_BIOM_PATHWAYS"].strip()
        else:
            # ... existing code ... 