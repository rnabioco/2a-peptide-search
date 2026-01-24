#!/usr/bin/env python3
"""
Generate sequence logos using Skylign API.

Accepts Stockholm alignment or HMM file and creates PNG logo via Skylign web service.
"""

import time
import gzip
import click
import requests
from pathlib import Path


def submit_to_skylign(file_path, api_url="http://skylign.org", processing="hmm", retry_attempts=3, retry_delay=5):
    """
    Submit alignment or HMM to Skylign API and retrieve logo UUID.

    Args:
        file_path: Path to Stockholm alignment or HMM file
        api_url: Base URL for Skylign API
        processing: Processing mode ('hmm' or 'obs')
        retry_attempts: Number of retry attempts on failure
        retry_delay: Delay in seconds between retries

    Returns:
        UUID string for retrieving the logo
    """
    url = api_url

    # Determine file format from extension
    file_path = Path(file_path)

    # Check if gzipped
    is_gzipped = file_path.suffix == '.gz'
    if is_gzipped:
        # Get the actual format from the extension before .gz
        actual_suffix = file_path.suffixes[-2] if len(file_path.suffixes) >= 2 else ''
    else:
        actual_suffix = file_path.suffix

    if actual_suffix == '.hmm':
        format_type = 'hmm'
    elif actual_suffix in ['.sto', '.stockholm']:
        format_type = 'stockholm'
    else:
        raise ValueError(f"Unsupported file format: {actual_suffix}")

    # Read file content (handle gzipped files)
    if is_gzipped:
        with gzip.open(file_path, 'rt') as f:
            file_content = f.read()
    else:
        with open(file_path, 'r') as f:
            file_content = f.read()

    # Prepare request
    data = {
        'processing': processing,
    }
    files = {
        'file': (file_path.name, file_content)
    }
    headers = {
        'Accept': 'application/json'
    }

    # Submit with retries
    for attempt in range(retry_attempts):
        try:
            response = requests.post(url, data=data, files=files, headers=headers, timeout=30)
            response.raise_for_status()

            # Parse UUID from response
            result = response.json()
            if 'url' in result:
                # Extract UUID from URL (format: /logo/UUID)
                uuid = result['url'].split('/')[-1]
                return uuid
            elif 'uuid' in result:
                return result['uuid']
            else:
                raise ValueError(f"No UUID in response: {result}")

        except requests.exceptions.RequestException as e:
            if attempt < retry_attempts - 1:
                click.echo(f"Attempt {attempt + 1} failed: {e}. Retrying in {retry_delay}s...")
                time.sleep(retry_delay)
            else:
                raise click.ClickException(f"Failed to submit to Skylign after {retry_attempts} attempts: {e}")


def download_logo(uuid, output_path, api_url="http://skylign.org", retry_attempts=10, retry_delay=3):
    """
    Download PNG logo from Skylign using UUID.

    Args:
        uuid: UUID returned from submission
        output_path: Path to save PNG logo
        api_url: Base URL for Skylign API
        retry_attempts: Number of retry attempts (increased for async processing)
        retry_delay: Delay in seconds between retries
    """
    url = f"{api_url}/logo/{uuid}"
    headers = {
        'Accept': 'image/png'
    }

    # Initial delay to let Skylign process the HMM
    time.sleep(2)

    for attempt in range(retry_attempts):
        try:
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()

            # Check if response is actually a PNG (starts with PNG magic bytes)
            # PNG files start with: 0x89 0x50 0x4E 0x47 (‰PNG)
            if len(response.content) < 8 or response.content[:4] != b'\x89PNG':
                # Not a PNG - job might still be processing
                error_text = response.content.decode('utf-8', errors='ignore')
                if attempt < retry_attempts - 1:
                    click.echo(f"    Job not ready (attempt {attempt + 1}): {error_text[:50]}. Retrying in {retry_delay}s...")
                    time.sleep(retry_delay)
                    continue
                else:
                    raise ValueError(f"Skylign did not return a PNG: {error_text[:100]}")

            # Save PNG content
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_path, 'wb') as f:
                f.write(response.content)

            return

        except requests.exceptions.RequestException as e:
            if attempt < retry_attempts - 1:
                click.echo(f"    Download attempt {attempt + 1} failed: {e}. Retrying in {retry_delay}s...")
                time.sleep(retry_delay)
            else:
                raise click.ClickException(f"Failed to download logo after {retry_attempts} attempts: {e}")


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('output_png', type=click.Path())
@click.option('--api-url', default='http://skylign.org', help='Skylign API base URL')
@click.option('--processing', default='hmm', type=click.Choice(['hmm', 'obs']),
              help='Processing method (hmm or obs)')
@click.option('--retry-attempts', default=3, help='Number of retry attempts')
@click.option('--retry-delay', default=5, help='Delay between retries (seconds)')
def main(input_file, output_png, api_url, processing, retry_attempts, retry_delay):
    """
    Generate sequence logo from alignment or HMM using Skylign API.

    INPUT_FILE: Stockholm alignment (.sto) or HMM file (.hmm)
    OUTPUT_PNG: Path to save PNG logo
    """
    click.echo(f"Submitting {input_file} to Skylign...")

    # Submit to Skylign
    uuid = submit_to_skylign(
        input_file,
        api_url=api_url,
        processing=processing,
        retry_attempts=retry_attempts,
        retry_delay=retry_delay
    )

    click.echo(f"Received UUID: {uuid}")

    # Download logo
    click.echo(f"Downloading logo...")
    download_logo(
        uuid,
        output_png,
        api_url=api_url,
        retry_attempts=retry_attempts,
        retry_delay=retry_delay
    )

    click.echo(f"Logo saved to {output_png}")


if __name__ == '__main__':
    main()
