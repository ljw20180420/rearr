<script setup>
import { ref } from 'vue';
import axios from 'axios';
const props = defineProps({
    params: {
        type: Array,
        default: []
    },
    descriptions: {
        type: Array,
        default: []
    },
    values: {
        type: Array,
        default: []
    },
    selectParams: {
        type: Array,
        default: []
    },
    selectDescriptions: {
        type: Array,
        default: []
    },
    selectValues: {
        type: Array,
        default: []
    }
});
const values = ref(props.values);
const selected = ref([]);
for (let i = 0; i < props.selectParams.length; ++i) {
    selected.value[i] = props.selectValues[i][0];
}
const status = ref("");
async function setValues() {
    try{
        const form = {};
        for (let i = 0; i < props.params.length; ++i) {
            form[props.params[i]] = values.value[i];
        }
        for (let i = 0; i < props.selectParams.length; ++i) {
            form[props.selectParams[i]] = selected.value[i];
        }
        status.value = "setting";
        const response = await axios.putForm("/setValues", form);
        status.value = response.data;
    } catch(error) {
        alert(error);
    }
}
function notAllSet() {
    for (let i = 0; i < values.value.length; ++i) {
        if (values.value[i] == null) {
            return true;
        }
    }
    return false;
}
</script>

<template>
    <div v-for="(param, i) in params">
        <span>{{ param }}</span>
        <input v-model="values[i]">
    </div>
    <div v-for="(selectDescription, i) in selectDescriptions">
        <span>{{ selectDescription }}</span>
        <select v-model="selected[i]">
            <option v-for="selectValue in selectValues[i]">
                {{ selectValue }}
            </option>
        </select>
    </div>
    <button @click="setValues" :disabled="notAllSet()">set values</button>
    <span>{{ status }}</span>
</template>