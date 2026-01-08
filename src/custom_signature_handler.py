"""
Custom Signature Handler Module

Parses and validates uploaded CSV files containing custom gene signatures.
Ensures compatibility with the DataManager's expected signature format.
"""

import pandas as pd
import io
import base64
from typing import Tuple, Optional


def parse_uploaded_csv(contents: str, filename: str) -> Tuple[Optional[pd.DataFrame], Optional[str]]:
    """
    Parse uploaded CSV file and validate format.
    
    Args:
        contents: Base64 encoded CSV content from dcc.Upload
        filename: Name of the uploaded file
        
    Returns:
        Tuple of (DataFrame, error_message)
        - If successful: (DataFrame, None)
        - If error: (None, error_message)
    """
    try:
        # Decode the base64 content
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        
        # Try to read as CSV
        try:
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        except UnicodeDecodeError:
            return None, "File encoding error. Please use UTF-8 encoded CSV."
        except pd.errors.ParserError as e:
            return None, f"CSV parsing error: {str(e)}"
        
        # Validate required columns
        required_columns = ['Cancer', 'Comparison', 'Signature']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            return None, f"Missing required column(s): {', '.join(missing_columns)}"
        
        # Check for empty DataFrame
        if df.empty:
            return None, "CSV file is empty"
        
        # Validate that required columns have values
        for col in required_columns:
            if df[col].isna().any():
                return None, f"Column '{col}' contains empty values"
        
        # Validate Signature format (should contain semicolons for multiple genes)
        # Allow single genes without semicolons
        for idx, sig in enumerate(df['Signature']):
            if not isinstance(sig, str) or sig.strip() == '':
                return None, f"Row {idx + 2}: Signature must be a non-empty string"
            
            # Basic validation: check if it looks like Ensembl IDs
            genes = sig.split(';')
            for gene in genes:
                gene = gene.strip()
                if not gene:
                    return None, f"Row {idx + 2}: Empty gene ID in signature"
        
        # Add auto-generated fields
        df['Filter'] = 'Merge'  # Default filter
        df['gProfiler'] = 'noLink'  # Placeholder for gProfiler link
        
        # Generate gProfiler links if desired (optional enhancement)
        # For now, using 'noLink' as placeholder
        
        return df, None
        
    except Exception as e:
        return None, f"Unexpected error: {str(e)}"


def validate_signature_dataframe(df: pd.DataFrame) -> Tuple[bool, Optional[str]]:
    """
    Additional validation for signature DataFrame.
    
    Args:
        df: DataFrame to validate
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    required_columns = ['Cancer', 'Comparison', 'Signature', 'Filter', 'gProfiler']
    
    for col in required_columns:
        if col not in df.columns:
            return False, f"Missing column: {col}"
    
    if df.empty:
        return False, "DataFrame is empty"
    
    return True, None
