import { test, expect } from '@playwright/test';

test('App should load and display pipeline overview', async ({ page }) => {
  await page.goto('/');
  await expect(page.locator('text=Pipeline Overview')).toBeVisible();
});
