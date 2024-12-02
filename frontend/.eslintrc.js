module.exports = {
    root: true,
    env: {
        browser: true,
        es2021: true,
        node: true,
    },
    extends: ['plugin:vue/recommended', 'eslint:recommended', '@vue/typescript/recommended', 'prettier'],
    ignorePatterns: ['**/*.min.js'],
    parserOptions: {
        ecmaVersion: 2020,
    },
    rules: {
        'vue/no-v-html': 'off',
        // TODO get rid of these rules for better code style
        'vue/require-default-prop': 'off',
        'vue/no-v-text-v-html-on-component': 'off',
        'vue/no-reserved-component-names': 'off',
        'vue/multi-word-component-names': 'off',
        '@typescript-eslint/no-explicit-any': 'off',
        '@typescript-eslint/no-inferrable-types': 'off',
        '@typescript-eslint/no-empty-interface': 'off',
    },
};
