<script setup>
const props = defineProps({
    params: {
        type: Array,
        default: []
    },
    selectParams: {
        type: Array,
        default: []
    },
    selectValues: {
        type: Array,
        default: []
    }
});
const pageStatus = defineModel();
</script>

<template>
    <div v-for="param, i in params">
        <span>{{ param }}</span>
        <input v-model="pageStatus[param]">
    </div>
    <div v-for="(selectParam, i) in selectParams">
        <span>{{ selectParam }}</span>
        <select v-model="pageStatus[selectParam]">
            <option v-for="selectValue in selectValues[i]">
                {{ selectValue }}
            </option>
        </select>
    </div>
</template>